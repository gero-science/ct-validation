#!/usr/bin/env python3
"""Aggregate clinical trials data from standardized parser outputs.

Discovers sources in source_dir/<source>_gene_drug.parquet (+ <source>_drug_indication.parquet).
Uses Union-Find drug clustering and per-source join keys.
Uses DuckDB for fast aggregation and joins.

Creates:
1. drug_id_mapping.parquet - Unified drug ID mapping
2. gene_drug_mapping.parquet - Combined gene-drug from all sources
3. drug_indication_mapping.parquet - Combined drug-indication
4. gene_indication_max_phase.parquet - Gene-indication with max phase
5. gene_drug_by_source.parquet - As (2), one row per source instead of collapsed
6. drug_indication_by_source.parquet - As (3), one row per source instead of collapsed

Filtering a combined table by source name is wrong: collapsing takes MAX(phase) and
`&`-joins source, so other sources' phases and action types land in that arm.
"""

import argparse
import logging
from pathlib import Path

import duckdb
import pandas as pd
import trial_status
import yaml

log = logging.getLogger(__name__)

# --- Per-source config ---

SOURCE_JOIN_KEYS = {
    "chembl": "chembl_id",
    "trialpanorama": "drug_name_normalized",
    "stitch": "stitch_id",
    "dgidb": "drug_name_normalized",
    "opentargets": "chembl_id",
}

# ChEMBL indications carry no trial status, and `phase_concluded` is a MAX across
# sources: one status-less record at a pair's top phase suppresses censoring for the
# whole pair. Only the indications are dropped — these drugs still feed the clustering.
DRUG_INDICATION_EXCLUDED_SOURCES = {"chembl"}

SOURCE_DISPLAY_NAMES = {
    "chembl": "ChEMBL",
    "trialpanorama": "TrialPanorama",
    "stitch": "STITCH",
    "dgidb": "DGIdb",
    "opentargets": "OpenTargets",
}

SCORE_COLUMNS = [
    "pchembl",
    "stitch_score",
    "dgidb_interaction_score",
    "dgidb_drug_specificity_score",
    "dgidb_gene_specificity_score",
    "dgidb_evidence_score",
]

ID_COLUMNS = ["chembl_id", "drugbank_id", "drug_name_normalized", "stitch_id"]
DRUG_METADATA_COLUMNS = ["molecule_type", "target_class", "trade_names", "synonyms"]

# Whose name names a drug cluster: curated names beat registry text, which is often
# a trial arm ("1.25 mg intravitreal bevacizumab") rather than a drug.
SOURCE_NAME_PRIORITY = ["ChEMBL", "OpenTargets", "DrugBank", "DGIdb", "TrialPanorama", "STITCH"]


def _drug_name_key(alias: str = "") -> str:
    """SQL for the normalized drug name, NULL when blank.

    Cluster membership is built with this expression and every join recomputes it,
    so any drift between the two silently drops rows.
    """
    col = f"{alias}.drug_name" if alias else "drug_name"
    return f"CASE WHEN TRIM(LOWER(COALESCE({col}, ''))) = '' THEN NULL ELSE TRIM(LOWER({col})) END"


def _join_key_sql(alias: str, join_key: str) -> tuple[str, str]:
    """Return (key expression, non-empty check) for joining `alias` to a cluster."""
    if join_key == "drug_name_normalized":
        return _drug_name_key(alias), f"{_drug_name_key(alias)} IS NOT NULL"
    return f"{alias}.{join_key}", f"{alias}.{join_key} IS NOT NULL AND {alias}.{join_key} != ''"


def load_dgidb_synonyms(drugs_path: Path) -> dict[str, str]:
    """Load DGIdb drugs.tsv and build a synonym -> canonical name mapping.

    Uses "Primary Drug Name" rows where drug_claim_name differs from drug_name,
    giving us trade names and alternative spellings mapped to canonical forms.
    All names are lowercased for case-insensitive matching.

    Returns:
        dict mapping lowercased synonym -> lowercased canonical name.
    """
    df = pd.read_csv(drugs_path, sep="\t", comment="#", low_memory=False)
    df = df[df["nomenclature"] == "Primary Drug Name"]

    synonym = df["drug_claim_name"].str.lower().str.strip()
    canonical = df["drug_name"].str.lower().str.strip()

    # Only keep where synonym differs from canonical (true synonyms)
    mask = synonym != canonical
    mapping = dict(zip(synonym[mask], canonical[mask]))
    log.info(f"Loaded {len(mapping):,} DGIdb drug synonyms")
    return mapping


class _UnionFind:
    """Lightweight Union-Find with path compression."""

    def __init__(self):
        self._parent = {}

    def find(self, x):
        if x not in self._parent:
            self._parent[x] = x
        if self._parent[x] != x:
            self._parent[x] = self.find(self._parent[x])
        return self._parent[x]

    def union(self, x, y):
        self._parent[self.find(x)] = self.find(y)

    def union_all(self, ids: list):
        """Union all identifiers together."""
        if len(ids) > 1:
            for j in range(1, len(ids)):
                self.union(ids[0], ids[j])
        elif len(ids) == 1:
            self.find(ids[0])

    def cluster_ids(self) -> dict:
        """Return mapping of root -> sequential cluster ID."""
        id_to_cluster = {}
        cluster_id = 0
        for identifier in self._parent:
            root = self.find(identifier)
            if root not in id_to_cluster:
                id_to_cluster[root] = cluster_id
                cluster_id += 1
        return id_to_cluster


def _collect_row_ids(
    id_arrays: dict[str, object],
    row: int,
    drug_synonyms: dict[str, str],
) -> list:
    """Collect all identifiers for a drug record row, including synonym links."""
    ids = [id_arrays[col][row] for col in ID_COLUMNS if pd.notna(id_arrays[col][row])]
    name = id_arrays["drug_name_normalized"][row]
    if pd.notna(name) and name in drug_synonyms:
        canonical = drug_synonyms[name]
        if canonical not in ids:
            ids.append(canonical)
    return ids


def create_drug_clusters(
    df: pd.DataFrame,
    drug_synonyms: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Assign cluster IDs to drugs using Union-Find algorithm.

    Args:
        df: Drug records with ID_COLUMNS.
        drug_synonyms: Optional mapping of synonym -> canonical name (lowercased).
            When provided, drugs whose drug_name_normalized matches a synonym are
            also unioned with the canonical name, bridging clusters via exact
            name matches.
    """
    uf = _UnionFind()
    synonyms = drug_synonyms or {}
    id_arrays = {col: df[col].to_numpy() for col in ID_COLUMNS}
    n_rows = len(df)

    for i in range(n_rows):
        uf.union_all(_collect_row_ids(id_arrays, i, synonyms))

    id_to_cluster = uf.cluster_ids()
    next_cluster = max(id_to_cluster.values(), default=-1) + 1

    clusters = []
    for i in range(n_rows):
        ids = _collect_row_ids(id_arrays, i, synonyms)
        if ids:
            clusters.append(id_to_cluster[uf.find(ids[0])])
        else:
            clusters.append(next_cluster)
            next_cluster += 1

    df["drug_cluster"] = clusters
    return df


def _get_connection() -> duckdb.DuckDBPyConnection:
    """Create a DuckDB connection with optimal settings."""
    return duckdb.connect(":memory:")


def _load_sources(
    source_dir: Path,
) -> tuple[dict[str, Path], dict[str, Path]]:
    """Discover gene_drug / drug_indication parquet paths per source."""
    log.info("Discovering source data...")
    gene_drug_paths = {}
    drug_indication_paths = {}

    for source_name in SOURCE_JOIN_KEYS:
        gd_path = source_dir / f"{source_name}_gene_drug.parquet"
        if gd_path.exists():
            log.info(f"  {source_name}_gene_drug: found")
            gene_drug_paths[source_name] = gd_path
        else:
            log.warning(f"  {source_name}_gene_drug: NOT FOUND")

        di_path = source_dir / f"{source_name}_drug_indication.parquet"
        if di_path.exists():
            log.info(f"  {source_name}_drug_indication: found")
            drug_indication_paths[source_name] = di_path

    return gene_drug_paths, drug_indication_paths


def _build_drug_id_mapping(
    conn: duckdb.DuckDBPyConnection,
    gene_drug_paths: dict[str, Path],
    drug_indication_paths: dict[str, Path],
    drug_synonyms: dict[str, str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Collect drug records from all sources and build unified drug ID mapping.

    Returns:
        (drug_id_mapping, membership): one display row per cluster, plus the
        identifier -> cluster lookup sources join through (see _cluster_membership).
    """
    log.info("Creating drug ID mapping...")

    # Build UNION ALL query to collect drug records from all sources.
    # Include both gene_drug and drug_indication files — a source like opentargets
    # may have drugs in one that aren't in the other.
    all_source_files = [
        *((name, path) for name, path in gene_drug_paths.items()),
        *((name, path) for name, path in drug_indication_paths.items()),
    ]
    # Deduplicate by path (same file shouldn't be scanned twice)
    seen_paths: set[Path] = set()
    unique_source_files = []
    for name, path in all_source_files:
        if path not in seen_paths:
            seen_paths.add(path)
            unique_source_files.append((name, path))
    union_parts = []
    for source_name, path in unique_source_files:
        display_name = SOURCE_DISPLAY_NAMES.get(source_name, source_name)
        # Read schema to check available columns
        schema_df = conn.execute(f"SELECT * FROM '{path}' LIMIT 0").df()
        cols = schema_df.columns.tolist()

        optional = ["chembl_id", "drugbank_id", "stitch_id", *DRUG_METADATA_COLUMNS]
        select_cols = ["drug_name", *(col if col in cols else f"NULL AS {col}" for col in optional)]

        union_parts.append(
            f"SELECT {', '.join(select_cols)}, '{display_name}' AS source FROM '{path}'",
        )

    if not union_parts:
        return pd.DataFrame(), pd.DataFrame()

    union_query = " UNION ALL ".join(union_parts)

    # Add normalized drug name and collect all records
    rank_cases = " ".join(f"WHEN '{name}' THEN {i}" for i, name in enumerate(SOURCE_NAME_PRIORITY))
    all_records = conn.execute(
        # sql
        f"""
        SELECT *,
            {_drug_name_key()} AS drug_name_normalized,
            CASE source {rank_cases} ELSE {len(SOURCE_NAME_PRIORITY)} END AS source_rank
        FROM ({union_query})
    """,
    ).df()

    log.info(f"  Total drug records: {len(all_records):,}")

    # Ensure all ID columns exist
    for col in ID_COLUMNS:
        if col not in all_records.columns:
            all_records[col] = None
    for col in DRUG_METADATA_COLUMNS:
        if col not in all_records.columns:
            all_records[col] = None

    # Run Union-Find clustering (Python)
    all_records = create_drug_clusters(all_records, drug_synonyms=drug_synonyms)
    log.info(f"  Drug clusters: {all_records['drug_cluster'].nunique():,}")

    # Register and aggregate with DuckDB
    conn.register("drug_records", all_records)

    # Unordered FIRST picks an arbitrary row of a parallel scan, which made this table
    # and every count downstream of it differ between runs. Name and normalized name
    # share one ORDER BY so they can never describe two different drugs.
    name_order = "ORDER BY source_rank, drug_name"
    result = conn.execute(
        # sql
        f"""
        SELECT
            drug_cluster AS our_drug_id,
            FIRST(drug_name {name_order}) FILTER (WHERE drug_name IS NOT NULL) AS drug_name,
            FIRST(drug_name_normalized {name_order})
                FILTER (WHERE drug_name IS NOT NULL) AS drug_name_normalized,
            MIN(chembl_id) AS chembl_id,
            MIN(drugbank_id) AS drugbank_id,
            MIN(stitch_id) AS stitch_id,
            MIN(molecule_type) AS molecule_type,
            MIN(target_class) AS target_class,
            MIN(trade_names) AS trade_names,
            MIN(synonyms) AS synonyms,
            STRING_AGG(DISTINCT source, '&' ORDER BY source) AS source
        FROM drug_records
        GROUP BY drug_cluster
        ORDER BY drug_cluster
    """,
    ).df()

    conn.unregister("drug_records")
    return result, _cluster_membership(all_records)


def _cluster_membership(all_records: pd.DataFrame) -> pd.DataFrame:
    """Map every identifier in a cluster to its cluster, as (id_type, id_value).

    Sources must join through this, not through drug_id_mapping: that table keeps one
    representative identifier per cluster, so a row keyed on any other member
    identifier silently fails to join and is dropped.
    """
    parts = [
        all_records[[col, "drug_cluster"]]
        .dropna(subset=[col])
        .rename(columns={col: "id_value", "drug_cluster": "our_drug_id"})
        .assign(id_type=col)
        for col in ID_COLUMNS
    ]
    membership = pd.concat(parts, ignore_index=True)
    # dropna() misses empty identifiers, which every record carrying one shares — they
    # would all collapse into a single cluster.
    membership = membership[membership["id_value"].astype(str).str.strip() != ""]
    membership = membership.drop_duplicates()

    # _UnionFind keys on the bare value with no id_type namespace, so a value used as
    # two kinds of identifier would merge two unrelated drugs.
    ambiguous = membership.groupby("id_value")["id_type"].nunique()
    if (clashes := int((ambiguous > 1).sum())) > 0:
        raise AssertionError(f"{clashes} identifiers are used as more than one id_type")

    log.info(f"  Cluster identifiers: {len(membership):,}")
    return membership


def _build_gene_drug(
    conn: duckdb.DuckDBPyConnection,
    gene_drug_paths: dict[str, Path],
    membership: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge each gene-drug source with drug clusters, concat and aggregate.

    Returns the collapsed table and the same edges kept one row per source.
    """
    log.info("Creating gene-drug mapping...")

    conn.register("drug_membership", membership)

    # Build UNION ALL with source-specific join logic
    union_parts = []
    for source_name, path in gene_drug_paths.items():
        display_name = SOURCE_DISPLAY_NAMES.get(source_name, source_name)
        join_key = SOURCE_JOIN_KEYS[source_name]

        # Read schema
        schema_df = conn.execute(f"SELECT * FROM '{path}' LIMIT 0").df()
        cols = schema_df.columns.tolist()

        # Build select with available columns
        optional_cols = [
            "subsource",
            "action_type",
            "drug_name",
            "chembl_id",
            "drugbank_id",
            "stitch_id",
            "molecule_type",
            *SCORE_COLUMNS,
        ]

        select_parts = ["gd.gene", "dm.our_drug_id", f"'{display_name}' AS source"]
        for col in optional_cols:
            if col in cols:
                select_parts.append(f"gd.{col}")
            else:
                select_parts.append(f"NULL AS {col}")

        key_expr, null_check = _join_key_sql("gd", join_key)

        union_parts.append(
            # sql
            f"""
            SELECT {", ".join(select_parts)}
            FROM '{path}' gd
            INNER JOIN drug_membership dm
                ON dm.id_type = '{join_key}' AND dm.id_value = {key_expr}
            WHERE gd.gene IS NOT NULL AND TRIM(gd.gene) != '' AND TRIM(gd.gene) != 'UNKNOWN'
                AND {null_check}
        """,
        )

    if not union_parts:
        conn.unregister("drug_membership")
        return pd.DataFrame(), pd.DataFrame()

    union_query = " UNION ALL ".join(union_parts)
    conn.execute(f"CREATE OR REPLACE TEMP TABLE gd_records AS {union_query}")

    # Build aggregation - check which score columns are available
    score_aggs = [f"MAX({col}) AS {col}" for col in SCORE_COLUMNS]

    result = conn.execute(
        # sql
        f"""
        SELECT
            gene,
            our_drug_id,
            STRING_AGG(DISTINCT source, '&' ORDER BY source) AS source,
            STRING_AGG(DISTINCT subsource, '&' ORDER BY subsource)
                FILTER (WHERE subsource IS NOT NULL) AS subsource,
            -- Collapsing to one value resolves every activation/inhibition conflict
            -- in the same direction; keep them all instead.
            STRING_AGG(DISTINCT action_type, '&' ORDER BY action_type)
                FILTER (WHERE action_type IS NOT NULL) AS action_type,
            MIN(drug_name) AS drug_name,
            MIN(chembl_id) AS chembl_id,
            MIN(drugbank_id) AS drugbank_id,
            MIN(stitch_id) AS stitch_id,
            MIN(molecule_type) AS molecule_type,
            {", ".join(score_aggs)}
        FROM gd_records
        GROUP BY gene, our_drug_id
    """,
    ).df()

    by_source = conn.execute(
        # sql
        """
        SELECT
            gene,
            our_drug_id,
            source,
            STRING_AGG(DISTINCT subsource, '&' ORDER BY subsource)
                FILTER (WHERE subsource IS NOT NULL) AS subsource,
            STRING_AGG(DISTINCT action_type, '&' ORDER BY action_type)
                FILTER (WHERE action_type IS NOT NULL) AS action_type
        FROM gd_records
        GROUP BY gene, our_drug_id, source
    """,
    ).df()

    conn.execute("DROP TABLE gd_records")
    conn.unregister("drug_membership")
    log.info(f"  Gene-drug pairs: {len(result):,}, genes: {result['gene'].nunique():,}")
    log.info(f"  Gene-drug pairs by source: {len(by_source):,}")
    return result, by_source


def _build_drug_indication(
    conn: duckdb.DuckDBPyConnection,
    drug_indication_paths: dict[str, Path],
    membership: pd.DataFrame,
    *,
    unknown_ongoing: bool = False,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Merge each drug-indication source with drug clusters, concat and aggregate.

    unknown_ongoing reclassifies Unknown-status records as still running. It is the
    sensitivity arm for reading them as concluded, not an alternative default.

    Returns the collapsed table and the same records kept one row per source.
    """
    log.info("Creating drug-indication mapping...")

    if not drug_indication_paths:
        log.warning("  No drug-indication data found")
        return pd.DataFrame(), pd.DataFrame()

    conn.register("drug_membership", membership)

    union_parts = []
    for source_name, path in drug_indication_paths.items():
        display_name = SOURCE_DISPLAY_NAMES.get(source_name, source_name)
        join_key = SOURCE_JOIN_KEYS[source_name]

        # Read schema
        schema_df = conn.execute(f"SELECT * FROM '{path}' LIMIT 0").df()
        cols = schema_df.columns.tolist()

        # Build select
        select_parts = ["dm.our_drug_id", f"'{display_name}' AS source"]
        for col in [
            "mesh_id",
            "mesh_heading",
            "efo_id",
            "efo_term",
            "phase",
            "subsource",
            "study_ids",
        ]:
            if col in cols:
                select_parts.append(f"di.{col}")
            else:
                select_parts.append(f"NULL AS {col}")
        # NULL reads as concluded below, so a status-less source silently un-censors
        # every pair it touches — see DRUG_INDICATION_EXCLUDED_SOURCES.
        if "is_concluded" not in cols:
            select_parts.append("CAST(NULL AS BOOLEAN) AS is_concluded")
        elif unknown_ongoing and "status" in cols:
            unknown = ", ".join(f"'{s}'" for s in sorted(trial_status.UNKNOWN_STATUSES))
            select_parts.append(
                f"CASE WHEN di.status IN ({unknown}) THEN FALSE ELSE di.is_concluded END"
                " AS is_concluded",
            )
        else:
            select_parts.append("di.is_concluded")
        # Handle reserved keyword separately
        if "references" in cols:
            select_parts.append('di."references" AS "references"')
        else:
            select_parts.append('NULL AS "references"')

        key_expr, null_check = _join_key_sql("di", join_key)

        union_parts.append(f"""
            SELECT {", ".join(select_parts)}
            FROM '{path}' di
            INNER JOIN drug_membership dm
                ON dm.id_type = '{join_key}' AND dm.id_value = {key_expr}
            WHERE {null_check}
        """)

    union_query = " UNION ALL ".join(union_parts)
    conn.execute(f"CREATE OR REPLACE TEMP TABLE di_records AS {union_query}")

    # An unmapped indication is not a pair: grouping them collapses every EFO-less
    # indication a drug has into one row whose phase spans all of them.
    unmapped = conn.execute("SELECT COUNT(*) FROM di_records WHERE efo_id IS NULL").fetchone()[0]
    log.info(f"  Dropping {unmapped:,} records with no EFO mapping")

    result = conn.execute(
        # sql
        """
        SELECT
            our_drug_id,
            efo_id,
            -- One EFO term routinely covers several MeSH headings, so id and heading
            -- share an ORDER BY: minimised apart they name two different diseases.
            FIRST(mesh_id ORDER BY mesh_id) FILTER (WHERE mesh_id IS NOT NULL) AS mesh_id,
            MIN(efo_term) AS efo_term,
            FIRST(mesh_heading ORDER BY mesh_id)
                FILTER (WHERE mesh_id IS NOT NULL) AS mesh_heading,
            MAX(phase) AS phase,
            MAX(phase) FILTER (WHERE COALESCE(is_concluded, TRUE)) AS phase_concluded,
            STRING_AGG(DISTINCT source, '&' ORDER BY source) AS source,
            STRING_AGG(DISTINCT subsource, '&' ORDER BY subsource)
                FILTER (WHERE subsource IS NOT NULL) AS subsource,
            MIN("references") AS "references",
            STRING_AGG(DISTINCT study_ids, '|' ORDER BY study_ids)
                FILTER (WHERE study_ids IS NOT NULL) AS study_ids
        FROM di_records
        WHERE efo_id IS NOT NULL
        GROUP BY our_drug_id, efo_id
    """,
    ).df()

    by_source = conn.execute(
        # sql
        """
        SELECT
            our_drug_id,
            efo_id,
            source,
            MAX(phase) AS phase,
            MAX(phase) FILTER (WHERE COALESCE(is_concluded, TRUE)) AS phase_concluded,
            STRING_AGG(DISTINCT subsource, '&' ORDER BY subsource)
                FILTER (WHERE subsource IS NOT NULL) AS subsource
        FROM di_records
        WHERE efo_id IS NOT NULL
        GROUP BY our_drug_id, efo_id, source
    """,
    ).df()

    conn.execute("DROP TABLE di_records")
    conn.unregister("drug_membership")
    log.info(f"  Drug-indication: {len(result):,}")
    log.info(f"  Drug-indication by source: {len(by_source):,}")
    return result, by_source


def is_ongoing_sql(prefix: str = "") -> str:
    """The undetermined-outcome rule as SQL, over an aggregate of one pair's records.

    Nothing concluded at the highest phase reached leaves the outcome undetermined. Approval
    (phase 4) is terminal: reaching it settles the outcome whatever post-marketing trials are
    still running.

    Not a valid global row filter — for a transition p->q a pair with max_phase >= q already
    succeeded and must still be counted. Censor only where `is_ongoing AND max_phase < q`.
    """
    return (
        f"COALESCE(MAX({prefix}phase) < 4, TRUE) "
        f"AND COALESCE(MAX({prefix}phase_concluded) < MAX({prefix}phase), TRUE)"
    )


def _build_gene_indication_max(
    conn: duckdb.DuckDBPyConnection,
    gene_drug: pd.DataFrame,
    drug_indication: pd.DataFrame,
) -> pd.DataFrame:
    """Join gene-drug with drug-indication and compute max phase per gene-indication."""
    log.info("Creating gene-indication-max phase...")

    if drug_indication.empty:
        return pd.DataFrame()

    conn.register("gene_drug", gene_drug)
    conn.register("drug_indication", drug_indication)

    # For source aggregation, we need to handle the '&'-delimited values
    # DuckDB can split and reaggregate
    result = conn.execute(
        # sql
        f"""
        WITH merged AS (
            SELECT
                gd.gene,
                di.mesh_id,
                di.mesh_heading,
                di.efo_id,
                di.efo_term,
                di.phase,
                di.phase_concluded,
                gd.source AS source_gd,
                di.source AS source_di,
                gd.pchembl,
                gd.stitch_score
            FROM gene_drug gd
            INNER JOIN drug_indication di ON gd.our_drug_id = di.our_drug_id
        ),
        -- Explode and reaggregate sources to handle already-combined values
        gd_sources_exploded AS (
            SELECT gene, efo_id,
                   UNNEST(STRING_SPLIT(source_gd, '&')) AS src
            FROM merged
        ),
        di_sources_exploded AS (
            SELECT gene, efo_id,
                   UNNEST(STRING_SPLIT(source_di, '&')) AS src
            FROM merged
        )
        SELECT
            m.gene,
            m.efo_id,
            MIN(m.mesh_id) AS mesh_id,
            MIN(m.mesh_heading) AS mesh_heading,
            MIN(m.efo_term) AS efo_term,
            MAX(m.phase) AS max_phase,
            MAX(m.phase_concluded) AS max_phase_concluded,
            {is_ongoing_sql("m.")} AS is_ongoing,
            (SELECT STRING_AGG(DISTINCT TRIM(src), '&' ORDER BY TRIM(src))
             FROM gd_sources_exploded g
             WHERE g.gene = m.gene
               AND g.efo_id IS NOT DISTINCT FROM m.efo_id) AS gene_drug_sources,
            (SELECT STRING_AGG(DISTINCT TRIM(src), '&' ORDER BY TRIM(src))
             FROM di_sources_exploded d
             WHERE d.gene = m.gene
               AND d.efo_id IS NOT DISTINCT FROM m.efo_id) AS drug_indication_sources,
            MAX(m.pchembl) AS pchembl,
            MAX(m.stitch_score) AS stitch_score
        FROM merged m
        GROUP BY m.gene, m.efo_id
    """,
    ).df()

    conn.unregister("gene_drug")
    conn.unregister("drug_indication")

    if result.empty:
        log.warning("  Gene-indication: 0 — gene-drug and drug-indication share no drugs")
        return result

    ongoing = int(result["is_ongoing"].sum())
    log.info(
        f"  Gene-indication: {len(result):,}, genes: {result['gene'].nunique():,}, "
        f"undetermined at highest phase: {ongoing:,} ({100 * ongoing / len(result):.1f}%)",
    )
    return result


def aggregate_clinical_trials(
    source_dir: Path,
    output_dir: Path,
    dgidb_drugs_path: Path | None = None,
    *,
    unknown_ongoing: bool = False,
):
    """Run full aggregation pipeline."""
    conn = _get_connection()

    # Load DGIdb drug synonyms if available
    drug_synonyms = None
    if dgidb_drugs_path and dgidb_drugs_path.exists():
        drug_synonyms = load_dgidb_synonyms(dgidb_drugs_path)

    gd_paths, di_paths = _load_sources(source_dir)
    if not gd_paths:
        log.error("No gene_drug sources found. Nothing to aggregate.")
        return

    drug_id_mapping, membership = _build_drug_id_mapping(conn, gd_paths, di_paths, drug_synonyms)
    gene_drug, gene_drug_by_source = _build_gene_drug(conn, gd_paths, membership)
    evidence_paths = {
        name: path
        for name, path in di_paths.items()
        if name not in DRUG_INDICATION_EXCLUDED_SOURCES
    }
    if dropped := sorted(di_paths.keys() & DRUG_INDICATION_EXCLUDED_SOURCES):
        log.info(f"Excluding from drug-indication evidence: {dropped}")
    if unknown_ongoing:
        log.info("Sensitivity arm: reading Unknown status as still running")
    drug_indication, drug_indication_by_source = _build_drug_indication(
        conn,
        evidence_paths,
        membership,
        unknown_ongoing=unknown_ongoing,
    )
    gene_indication_max = _build_gene_indication_max(conn, gene_drug, drug_indication)

    conn.close()

    output_dir.mkdir(parents=True, exist_ok=True)
    drug_id_mapping.to_parquet(output_dir / "drug_id_mapping.parquet", index=False)
    gene_drug.to_parquet(output_dir / "gene_drug_mapping.parquet", index=False)
    drug_indication.to_parquet(output_dir / "drug_indication_mapping.parquet", index=False)
    gene_indication_max.to_parquet(output_dir / "gene_indication_max_phase.parquet", index=False)
    gene_drug_by_source.to_parquet(output_dir / "gene_drug_by_source.parquet", index=False)
    drug_indication_by_source.to_parquet(
        output_dir / "drug_indication_by_source.parquet",
        index=False,
    )
    log.info(f"Saved to {output_dir}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--source-dir", type=Path, help="Directory containing source subdirs")
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Output directory (default: source_dir/aggregated)",
    )
    parser.add_argument("--dgidb-drugs", type=Path, help="DGIdb drugs.tsv for synonym matching")
    parser.add_argument(
        "--unknown-ongoing",
        action="store_true",
        help="Sensitivity arm: read Unknown trial status as still running, not concluded",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    config_path = (
        args.config or Path(__file__).parent.parent.parent.parent / "configs" / "parsing.yaml"
    )
    if config_path.exists():
        with config_path.open() as f:
            cfg = yaml.safe_load(f)
        ct_inputs = cfg.get("clinical_trials", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ct_inputs, output_dir = {}, Path()

    source_dir = args.source_dir or output_dir / "clinical_trials"
    # Its own directory: the sensitivity arm is not a drop-in replacement for the aggregate
    # the benchmark reads, and would otherwise silently overwrite it.
    default_name = "aggregated_unknown_ongoing" if args.unknown_ongoing else "aggregated"
    agg_output_dir = args.output_dir or source_dir / default_name
    dgidb_drugs = args.dgidb_drugs or (Path(p) if (p := ct_inputs.get("dgidb_drugs")) else None)

    aggregate_clinical_trials(
        source_dir,
        agg_output_dir,
        dgidb_drugs_path=dgidb_drugs,
        unknown_ongoing=args.unknown_ongoing,
    )


if __name__ == "__main__":
    main()
