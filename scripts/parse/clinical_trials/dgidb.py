#!/usr/bin/env python3
"""Parse DGIdb interactions to standardized gene-drug parquet.

Reads:
- DGIdb interactions TSV (interactions.tsv from https://dgidb.org/downloads)

Outputs:
- clinical_trials/dgidb_gene_drug.parquet
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def parse_dgidb(interactions_path: Path, output_path: Path) -> pd.DataFrame:
    """Parse DGIdb interactions to gene-drug mapping."""
    df = pd.read_csv(interactions_path, sep="\t", comment="#")
    log.info(f"Loaded interactions: {len(df):,} rows")

    # Extract chembl_id and drugbank_id from drug_concept_id
    concept = df["drug_concept_id"].fillna("")
    df["chembl_id"] = concept.where(concept.str.startswith("chembl:"), "").str.replace(
        "chembl:",
        "",
        regex=False,
    )
    df["drugbank_id"] = concept.where(concept.str.startswith("drugbank:"), "").str.replace(
        "drugbank:",
        "",
        regex=False,
    )
    df["chembl_id"] = df["chembl_id"].replace("", None)
    df["drugbank_id"] = df["drugbank_id"].replace("", None)

    result = df.rename(
        columns={
            "gene_name": "gene",
            "drug_concept_id": "dgidb_id",
            "interaction_types": "action_type",
            "interaction_score": "dgidb_interaction_score",
            "drug_specificity_score": "dgidb_drug_specificity_score",
            "gene_specificity_score": "dgidb_gene_specificity_score",
            "evidence_score": "dgidb_evidence_score",
        },
    )
    result["subsource"] = "DGIdb"

    keep = [
        "gene",
        "drug_name",
        "dgidb_id",
        "chembl_id",
        "drugbank_id",
        "subsource",
        "action_type",
        "dgidb_interaction_score",
        "dgidb_drug_specificity_score",
        "dgidb_gene_specificity_score",
        "dgidb_evidence_score",
    ]
    result = result[keep].drop_duplicates()

    log.info(
        f"Gene-drug: {len(result):,} rows, {result['gene'].nunique():,} genes, "
        f"{result['drug_name'].nunique():,} drugs",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--interactions", type=Path, help="DGIdb interactions TSV")
    parser.add_argument("--output", type=Path, help="Output parquet path")
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

    parse_dgidb(
        interactions_path=args.interactions or Path(ct_inputs.get("dgidb_interactions", "")),
        output_path=args.output or output_dir / "clinical_trials" / "dgidb_gene_drug.parquet",
    )


if __name__ == "__main__":
    main()
