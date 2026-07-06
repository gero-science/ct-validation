"""Match genetic evidence with clinical trials using semantic similarity."""

from pathlib import Path

import duckdb
import pandas as pd


def _setup_source(
    conn: duckdb.DuckDBPyConnection,
    name: str,
    data: pd.DataFrame | str | Path,
) -> None:
    """Register DataFrame or create view for parquet path."""
    if isinstance(data, pd.DataFrame):
        conn.register(name, data)
    else:
        conn.execute(f"CREATE OR REPLACE VIEW {name} AS SELECT * FROM read_parquet('{data}')")


def get_expanded_disease_set(
    efo_ids: set[str] | list[str],
    similarity_pairs: pd.DataFrame | str | Path,
    similarity_threshold: float,
) -> set[str]:
    """
    Expand disease IDs via semantic similarity.

    Returns set of efo_ids: input diseases + similar diseases (similarity >= threshold).
    If threshold >= 1.0, returns only input diseases (no expansion).

    Useful for filtering clinical trials to a disease neighborhood before validation.
    """
    direct_set = set(efo_ids)

    if similarity_threshold >= 1.0:
        return direct_set

    conn = duckdb.connect()
    _setup_source(conn, "sim", similarity_pairs)
    conn.register("input_diseases", pd.DataFrame({"efo_id": list(direct_set)}))

    # sql — single-direction probe. Correct because the lookup is required to be
    # symmetric (see README data contract); validate once with check_similarity_lookup.
    query = """
        SELECT DISTINCT sim.efo_id_2 AS efo_id
        FROM input_diseases d
        INNER JOIN sim ON d.efo_id = sim.efo_id_1 AND sim.similarity >= $threshold
        WHERE sim.efo_id_2 IS NOT NULL
    """
    expanded = conn.execute(query, {"threshold": similarity_threshold}).fetchall()
    conn.close()

    return direct_set | {row[0] for row in expanded}


def create_matched_pairs_set(
    genetic_evidence: pd.DataFrame | str | Path,
    clinical_trials: pd.DataFrame | str | Path,
    similarity_pairs: pd.DataFrame | str | Path | None,
    similarity_threshold: float,
    gene_universe: set[str] | None = None,
    *,
    con: duckdb.DuckDBPyConnection | None = None,
) -> set[tuple[str, str]]:
    """
    Find (gene, efo_id) pairs in clinical_trials that have matching genetic_evidence.

    Returns CT pairs where at least one GE match exists with similarity >= threshold.
    If similarity_pairs is None, performs exact matching on (gene, efo_id).
    If con is provided, it is reused (caller owns lifecycle).
    """
    owned = con is None
    conn = con or duckdb.connect()

    _setup_source(conn, "ge", genetic_evidence)
    _setup_source(conn, "ct", clinical_trials)

    if gene_universe:
        conn.register("gu", pd.DataFrame({"gene": list(gene_universe)}))
        gu_join = "INNER JOIN gu ON ct.gene = gu.gene"
    else:
        gu_join = ""

    if similarity_pairs is not None:
        _setup_source(conn, "sim", similarity_pairs)
        # sql — single-direction probe; correct on the required symmetric lookup
        # (validate once with check_similarity_lookup).
        query = f"""
            SELECT DISTINCT ct.gene, ct.efo_id
            FROM ct
            {gu_join}
            INNER JOIN sim ON ct.efo_id = sim.efo_id_2 AND sim.similarity >= $threshold
            INNER JOIN ge ON ct.gene = ge.gene AND sim.efo_id_1 = ge.efo_id
            WHERE ct.gene IS NOT NULL AND ct.efo_id IS NOT NULL
              AND ge.gene IS NOT NULL AND ge.efo_id IS NOT NULL
        """
        result = conn.execute(query, {"threshold": similarity_threshold}).fetchall()
    else:
        # sql
        query = f"""
            SELECT DISTINCT ct.gene, ct.efo_id
            FROM ct
            {gu_join}
            INNER JOIN ge ON ct.gene = ge.gene AND ct.efo_id = ge.efo_id
            WHERE ct.gene IS NOT NULL AND ct.efo_id IS NOT NULL
              AND ge.gene IS NOT NULL AND ge.efo_id IS NOT NULL
        """
        result = conn.execute(query).fetchall()

    if owned:
        conn.close()

    return set(result)


def create_matched_pairs_df(
    genetic_evidence: pd.DataFrame | str | Path,
    clinical_trials: pd.DataFrame | str | Path,
    similarity_pairs: pd.DataFrame | str | Path | None,
    similarity_threshold: float,
    gene_universe: set[str] | None = None,
    *,
    con: duckdb.DuckDBPyConnection | None = None,
) -> pd.DataFrame:
    """
    Return detailed target matches for audit and debugging.

    Columns: gene, ct_efo_id, ge_efo_id, similarity, match_type
    (match_type is ``exact`` when ct_efo_id equals ge_efo_id, else ``similarity``).
    """
    owned = con is None
    conn = con or duckdb.connect()

    _setup_source(conn, "ge", genetic_evidence)
    _setup_source(conn, "ct", clinical_trials)

    if gene_universe:
        conn.register("gu", pd.DataFrame({"gene": list(gene_universe)}))
        gu_join = "INNER JOIN gu ON ct.gene = gu.gene"
    else:
        gu_join = ""

    if similarity_pairs is not None:
        _setup_source(conn, "sim", similarity_pairs)
        # sql — single-direction probe; correct on the required symmetric lookup.
        query = f"""
            SELECT DISTINCT
                ct.gene,
                ct.efo_id AS ct_efo_id,
                ge.efo_id AS ge_efo_id,
                sim.similarity,
                CASE WHEN ct.efo_id = ge.efo_id THEN 'exact' ELSE 'similarity' END AS match_type
            FROM ct
            {gu_join}
            INNER JOIN sim ON ct.efo_id = sim.efo_id_2 AND sim.similarity >= $threshold
            INNER JOIN ge ON ct.gene = ge.gene AND sim.efo_id_1 = ge.efo_id
            WHERE ct.gene IS NOT NULL AND ct.efo_id IS NOT NULL
              AND ge.gene IS NOT NULL AND ge.efo_id IS NOT NULL
        """
        df = conn.execute(query, {"threshold": similarity_threshold}).df()
    else:
        query = f"""
            SELECT DISTINCT
                ct.gene,
                ct.efo_id AS ct_efo_id,
                ge.efo_id AS ge_efo_id,
                1.0 AS similarity,
                'exact' AS match_type
            FROM ct
            {gu_join}
            INNER JOIN ge ON ct.gene = ge.gene AND ct.efo_id = ge.efo_id
            WHERE ct.gene IS NOT NULL AND ct.efo_id IS NOT NULL
              AND ge.gene IS NOT NULL AND ge.efo_id IS NOT NULL
        """
        df = conn.execute(query).df()

    if owned:
        conn.close()

    return df.sort_values(["gene", "ct_efo_id", "ge_efo_id"]).reset_index(drop=True)


def check_similarity_lookup(
    similarity_pairs: pd.DataFrame | str | Path,
    *,
    require_diagonal: bool = True,
) -> None:
    """
    Assert a similarity lookup meets the symmetric-matrix contract; raise if not.

    Matching probes a single orientation of each pair (for speed), which is correct
    only when the table is symmetric: every off-diagonal pair ``(a, b)`` must also be
    stored as ``(b, a)``. When ``require_diagonal`` is set, every term must additionally
    have a self-similarity row ``(a, a)`` (exact matching relies on it). Run this once
    per lookup — after loading it, before a validation loop — so a malformed table fails
    loudly instead of letting matching silently drop half the disease neighbourhood.

    Raises:
        ValueError: if the table is not symmetric, or (when ``require_diagonal``) any
            term lacks a self-similarity row. The message lists a few offending entries.
    """
    conn = duckdb.connect()
    try:
        _setup_source(conn, "sim", similarity_pairs)

        # sql — off-diagonal pairs whose mirror (b, a) is absent
        asymmetric = conn.execute(
            """
            SELECT a.efo_id_1, a.efo_id_2
            FROM sim a
            ANTI JOIN sim b ON a.efo_id_1 = b.efo_id_2 AND a.efo_id_2 = b.efo_id_1
            WHERE a.efo_id_1 <> a.efo_id_2
            LIMIT 5
            """,
        ).fetchall()
        if asymmetric:
            examples = ", ".join(f"({a}, {b})" for a, b in asymmetric)
            raise ValueError(
                "similarity_lookup is not symmetric: some off-diagonal pairs are missing "
                f"their mirror (b, a), e.g. {examples}. Matching assumes a symmetric table; "
                "store both orientations (or regenerate via scripts/r/create_efo_similarity_matrix.R).",
            )

        if require_diagonal:
            # sql — terms that never appear in a self-similarity (a, a) row
            missing_diagonal = conn.execute(
                """
                SELECT t.efo_id
                FROM (SELECT efo_id_1 AS efo_id FROM sim UNION SELECT efo_id_2 FROM sim) t
                ANTI JOIN sim s ON s.efo_id_1 = t.efo_id AND s.efo_id_2 = t.efo_id
                LIMIT 5
                """,
            ).fetchall()
            if missing_diagonal:
                examples = ", ".join(row[0] for row in missing_diagonal)
                raise ValueError(
                    "similarity_lookup is missing diagonal (self-similarity) rows for some "
                    f"terms, which exact matching relies on, e.g. {examples}.",
                )
    finally:
        conn.close()
