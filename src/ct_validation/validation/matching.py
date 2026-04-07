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

    # sql
    query = """
        SELECT DISTINCT sim.efo_id_2
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
        # sql
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
