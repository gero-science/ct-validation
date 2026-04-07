#!/usr/bin/env python3
"""Parse pre-fetched ChEMBL data to standardized gene-drug and drug-indication parquets.

Reads:
- chembl_gene_drug: Pre-fetched gene-drug parquet (from scripts/fetch/chembl_fetch.py)
- chembl_indications: ChEMBL indications TSV (from ChEMBL website drug indications download)

Outputs:
- clinical_trials/chembl_gene_drug.parquet
- clinical_trials/chembl_drug_indication.parquet
"""

import argparse
import logging
from itertools import zip_longest
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def parse_gene_drug(gene_drug_path: Path) -> pd.DataFrame:
    """Parse pre-fetched gene-drug parquet to standard schema."""
    df = pd.read_parquet(gene_drug_path)
    log.info(f"Loaded gene-drug: {len(df):,} rows")

    result = df.rename(
        columns={
            "initial_target_name": "gene",
            "molecule_chembl_id": "chembl_id",
            "Parent Molecule Name": "drug_name",
            "Parent Molecule Type": "molecule_type",
        },
    )
    result["subsource"] = "ChEMBL"

    keep = ["gene", "drug_name", "chembl_id", "molecule_type", "subsource", "action_type"]
    keep = [c for c in keep if c in result.columns]
    result = result[keep].drop_duplicates()

    log.info(f"Gene-drug: {len(result):,} rows, {result['gene'].nunique():,} genes")
    return result


def parse_drug_indication(indications_path: Path) -> pd.DataFrame:
    """Parse ChEMBL indications TSV to standard schema with EFO explosion."""
    df = pd.read_csv(indications_path, sep="\t")
    log.info(f"Loaded indications: {len(df):,} rows")

    df = df.rename(
        columns={
            "MESH ID": "mesh_id",
            "MESH Heading": "mesh_heading",
            "Max Phase for Indication": "phase",
            "Parent Molecule ChEMBL ID": "chembl_id",
            "Parent Molecule Name": "drug_name",
            "Parent Molecule Type": "molecule_type",
            "References": "references",
        },
    )

    df["phase"] = pd.to_numeric(df["phase"], errors="coerce")
    df = df[(df["phase"].notna()) & (df["phase"] >= 0) & (df["phase"] <= 4)]
    df = df[df["phase"] == df["phase"].astype(int)]  # drop non-integer phases (e.g. 0.5)
    df["phase"] = df["phase"].astype("Int64")

    # Explode pipe-delimited EFO IDs/Terms (handles mismatched lengths)
    efo_ids_split = df["EFO IDs"].str.split("|")
    efo_terms_split = df["EFO Terms"].str.split("|")
    df["_efo_pairs"] = [
        list(zip_longest(ids or [None], terms or [None], fillvalue=None))
        for ids, terms in zip(
            efo_ids_split.where(efo_ids_split.notna(), [[None]]),
            efo_terms_split.where(efo_terms_split.notna(), [[None]]),
        )
    ]
    df = df.explode("_efo_pairs")
    df["efo_id"] = df["_efo_pairs"].str[0]
    df["efo_term"] = df["_efo_pairs"].str[1]

    df["subsource"] = "ChEMBL"

    keep = [
        "efo_id",
        "phase",
        "drug_name",
        "chembl_id",
        "mesh_id",
        "mesh_heading",
        "efo_term",
        "subsource",
        "references",
        "molecule_type",
    ]
    result = df[keep].drop_duplicates()

    log.info(f"Drug-indication: {len(result):,} rows, {result['chembl_id'].nunique():,} drugs")
    return result


def parse_chembl(
    gene_drug_path: Path,
    indications_path: Path | None,
    gene_drug_output: Path,
    drug_indication_output: Path | None,
) -> pd.DataFrame:
    """Parse ChEMBL data to standardized outputs."""
    gene_drug_output.parent.mkdir(parents=True, exist_ok=True)

    gene_drug = parse_gene_drug(gene_drug_path)
    gene_drug.to_parquet(gene_drug_output, index=False)
    log.info(f"Saved gene-drug to {gene_drug_output}")

    if indications_path and indications_path.exists():
        drug_indication = parse_drug_indication(indications_path)
        if drug_indication_output:
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
    parser.add_argument("--gene-drug", type=Path, help="Pre-fetched gene-drug parquet")
    parser.add_argument("--indications", type=Path, help="ChEMBL indications TSV")
    parser.add_argument("--output-dir", type=Path, help="Output directory")
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

    ct_out = output_dir / "clinical_trials"
    indications = args.indications or (
        Path(p) if (p := ct_inputs.get("chembl_indications")) else None
    )

    parse_chembl(
        gene_drug_path=args.gene_drug or Path(ct_inputs.get("chembl_gene_drug", "")),
        indications_path=indications,
        gene_drug_output=args.output_dir or ct_out / "chembl_gene_drug.parquet",
        drug_indication_output=ct_out / "chembl_drug_indication.parquet" if indications else None,
    )


if __name__ == "__main__":
    main()
