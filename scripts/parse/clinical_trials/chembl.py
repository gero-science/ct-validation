#!/usr/bin/env python3
"""Parse ChEMBL SQLite database to standardized gene-drug and drug-indication parquets.

Reads:
- ChEMBL SQLite database (downloaded via scripts/fetch/chembl_fetch.sh)

Outputs:
- clinical_trials/chembl_gene_drug.parquet
- clinical_trials/chembl_drug_indication.parquet
"""

import argparse
import logging
import sqlite3
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def _query_gene_drug(conn: sqlite3.Connection, min_pchembl: float) -> pd.DataFrame:
    """Extract gene-drug mappings: human single-protein targets with pChEMBL >= threshold."""
    log.info(f"Querying gene-drug mappings (pChEMBL >= {min_pchembl})...")
    # sql
    query = """
        SELECT DISTINCT
            cs.component_synonym AS gene,
            parent_md.pref_name AS drug_name,
            parent_md.chembl_id,
            parent_md.molecule_type
        FROM target_dictionary td
        JOIN target_components tc ON td.tid = tc.tid
        JOIN component_synonyms cs
            ON tc.component_id = cs.component_id AND cs.syn_type = 'GENE_SYMBOL'
        JOIN assays a ON td.tid = a.tid
        JOIN activities act ON a.assay_id = act.assay_id
        JOIN molecule_hierarchy mh ON act.molregno = mh.molregno
        JOIN molecule_dictionary parent_md ON mh.parent_molregno = parent_md.molregno
        WHERE td.target_type = 'SINGLE PROTEIN'
          AND td.organism = 'Homo sapiens'
          AND act.pchembl_value >= ?
          AND parent_md.max_phase > 0
    """
    df = pd.read_sql_query(query, conn, params=(min_pchembl,))
    df["subsource"] = "ChEMBL"
    log.info(f"Gene-drug: {len(df):,} rows, {df['gene'].nunique():,} genes")
    return df


def _query_drug_indication(conn: sqlite3.Connection) -> pd.DataFrame:
    """Extract drug-indication mappings with EFO terms."""
    log.info("Querying drug-indication mappings...")
    # sql
    query = """
        WITH refs AS (
            SELECT
                drugind_id,
                GROUP_CONCAT(
                    'Type: ' || ref_type || ' RefID: ' || ref_id || ' URL: ' || ref_url,
                    '|'
                ) AS "references"
            FROM indication_refs
            GROUP BY drugind_id
        )
        SELECT DISTINCT
            di.efo_id,
            CAST(di.max_phase_for_ind AS INTEGER) AS phase,
            parent_md.pref_name AS drug_name,
            parent_md.chembl_id,
            di.mesh_id,
            di.mesh_heading,
            di.efo_term,
            parent_md.molecule_type,
            r."references"
        FROM drug_indication di
        JOIN molecule_hierarchy mh ON di.molregno = mh.molregno
        JOIN molecule_dictionary parent_md ON mh.parent_molregno = parent_md.molregno
        LEFT JOIN refs r ON di.drugind_id = r.drugind_id
        WHERE di.max_phase_for_ind >= 0
          AND di.max_phase_for_ind = CAST(di.max_phase_for_ind AS INTEGER)
          AND di.efo_id IS NOT NULL
    """
    df = pd.read_sql_query(query, conn)
    # Convert EFO URIs to CURIEs:
    #   http://www.ebi.ac.uk/efo/EFO_0000249 -> EFO:0000249
    #   http://purl.obolibrary.org/obo/MONDO_0005015 -> MONDO:0005015
    df["efo_id"] = df["efo_id"].str.rsplit("/", n=1).str[-1].str.replace("_", ":", n=1)
    df["subsource"] = "ChEMBL"
    log.info(f"Drug-indication: {len(df):,} rows, {df['chembl_id'].nunique():,} drugs")
    return df


def parse_chembl(
    db_path: Path,
    gene_drug_output: Path,
    drug_indication_output: Path,
    min_pchembl: float = 7.0,
) -> pd.DataFrame:
    """Parse ChEMBL SQLite database to standardized outputs."""
    log.info(f"Opening ChEMBL database: {db_path}")
    conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)

    gene_drug = _query_gene_drug(conn, min_pchembl)
    drug_indication = _query_drug_indication(conn)
    conn.close()

    gene_drug_output.parent.mkdir(parents=True, exist_ok=True)
    gene_drug.to_parquet(gene_drug_output, index=False)
    log.info(f"Saved gene-drug to {gene_drug_output}")

    drug_indication_output.parent.mkdir(parents=True, exist_ok=True)
    drug_indication.to_parquet(drug_indication_output, index=False)
    log.info(f"Saved drug-indication to {drug_indication_output}")

    return gene_drug


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--db", type=Path, help="ChEMBL SQLite database path")
    parser.add_argument("--output-dir", type=Path, help="Output directory")
    parser.add_argument("--min-pchembl", type=float, help="Minimum pChEMBL threshold")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    config_path = (
        args.config
        or Path(__file__).parent.parent.parent.parent / "configs" / "parsing.yaml"
    )
    if config_path.exists():
        with config_path.open() as f:
            cfg = yaml.safe_load(f)
        ct_inputs = cfg.get("clinical_trials", {})
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ct_inputs, thresholds, output_dir = {}, {}, Path()

    out_dir = args.output_dir or output_dir / "clinical_trials"

    parse_chembl(
        db_path=args.db or Path(ct_inputs.get("chembl_db", "")),
        gene_drug_output=out_dir / "chembl_gene_drug.parquet",
        drug_indication_output=out_dir / "chembl_drug_indication.parquet",
        min_pchembl=args.min_pchembl or thresholds.get("chembl_min_pchembl", 7.0),
    )


if __name__ == "__main__":
    main()
