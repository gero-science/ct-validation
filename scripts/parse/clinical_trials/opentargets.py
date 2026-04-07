#!/usr/bin/env python3
"""Parse Open Targets known_drug data to standardized gene-drug and drug-indication parquets.

Reads:
- opentargets_known_drug: Open Targets known_drug parquet directory

Outputs:
- clinical_trials/opentargets_gene_drug.parquet
- clinical_trials/opentargets_drug_indication.parquet
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def _extract_efo_id(disease_id: str, ancestors: list | None) -> str | None:
    """Extract EFO ID from disease_id or ancestors list, converting _ to : format."""
    if pd.isna(disease_id):
        return None
    if disease_id.startswith("EFO_"):
        return disease_id.replace("EFO_", "EFO:")
    # Try to find EFO in ancestors
    if ancestors is not None:
        for anc in ancestors:
            if isinstance(anc, str) and anc.startswith("EFO_"):
                return anc.replace("EFO_", "EFO:")
    return None


def _join_list(val) -> str | None:
    """Join list/array to pipe-delimited string."""
    if val is None:
        return None
    if isinstance(val, float) and pd.isna(val):
        return None
    if hasattr(val, "__iter__") and not isinstance(val, str):
        items = [str(v) for v in val if v is not None and str(v)]
        return "|".join(items) if items else None
    return str(val) or None


def _extract_subsources(urls) -> str:
    """Extract source names from urls array."""
    if urls is None or (isinstance(urls, float) and pd.isna(urls)):
        return "OpenTargets"
    sources = set()
    try:
        for u in urls:
            if isinstance(u, dict) and "niceName" in u:
                sources.add(u["niceName"])
    except (TypeError, ValueError):
        pass
    return "|".join(sorted(sources)) if sources else "OpenTargets"


def parse_gene_drug(df: pd.DataFrame) -> pd.DataFrame:
    """Extract gene-drug pairs from Open Targets known_drug."""
    result = df[
        [
            "approvedSymbol",
            "prefName",
            "drugId",
            "drugType",
            "mechanismOfAction",
            "targetClass",
            "tradeNames",
            "synonyms",
            "urls",
        ]
    ].copy()

    # Convert list columns to pipe-delimited strings
    for col in ["targetClass", "tradeNames", "synonyms"]:
        result[col] = result[col].apply(_join_list)

    # Extract subsources from urls
    result["subsource"] = result["urls"].apply(_extract_subsources)

    result = result.rename(
        columns={
            "approvedSymbol": "gene",
            "prefName": "drug_name",
            "drugId": "chembl_id",
            "drugType": "molecule_type",
            "mechanismOfAction": "action_type",
            "targetClass": "target_class",
            "tradeNames": "trade_names",
        },
    )

    result = result[result["gene"].notna() & (result["gene"] != "")]
    result = result[result["chembl_id"].notna()]
    result = result.drop(columns=["urls"]).drop_duplicates()

    log.info(f"Gene-drug: {len(result):,} rows, {result['gene'].nunique():,} genes")
    return result


def parse_drug_indication(df: pd.DataFrame) -> pd.DataFrame:
    """Extract drug-indication pairs from Open Targets known_drug."""
    result = df[
        ["drugId", "prefName", "drugType", "diseaseId", "label", "ancestors", "phase", "urls"]
    ].copy()

    # Extract EFO IDs
    result["efo_id"] = result.apply(
        lambda r: _extract_efo_id(r["diseaseId"], r["ancestors"]),
        axis=1,
    )

    # Filter to rows with valid EFO and phase
    result["phase"] = pd.to_numeric(result["phase"], errors="coerce")
    result = result[(result["phase"].notna()) & (result["phase"] >= 1) & (result["phase"] <= 4)]
    result = result[result["phase"] == result["phase"].astype(int)]  # drop non-integer (0.5)
    result["phase"] = result["phase"].astype("Int64")

    # Extract subsources from urls
    result["subsource"] = result["urls"].apply(_extract_subsources)

    result = result.rename(
        columns={
            "drugId": "chembl_id",
            "prefName": "drug_name",
            "drugType": "molecule_type",
            "label": "efo_term",
        },
    )

    keep = ["efo_id", "phase", "drug_name", "chembl_id", "efo_term", "subsource", "molecule_type"]
    result = result[keep].drop_duplicates()

    log.info(
        f"Drug-indication: {len(result):,} rows, {result['chembl_id'].nunique():,} drugs, "
        f"{result['efo_id'].nunique():,} diseases",
    )
    return result


def parse_opentargets(
    known_drug_path: Path,
    gene_drug_output: Path,
    drug_indication_output: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse Open Targets known_drug to standardized outputs."""
    log.info(f"Loading Open Targets known_drug from {known_drug_path}")
    df = pd.read_parquet(known_drug_path)
    log.info(f"Loaded {len(df):,} rows")

    gene_drug_output.parent.mkdir(parents=True, exist_ok=True)
    drug_indication_output.parent.mkdir(parents=True, exist_ok=True)

    gene_drug = parse_gene_drug(df)
    gene_drug.to_parquet(gene_drug_output, index=False)

    drug_indication = parse_drug_indication(df)
    drug_indication.to_parquet(drug_indication_output, index=False)

    log.info(f"Saved to {gene_drug_output} and {drug_indication_output}")
    return gene_drug, drug_indication


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--known-drug", type=Path, help="Open Targets known_drug parquet dir")
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

    ct_out = args.output_dir or output_dir / "clinical_trials"
    parse_opentargets(
        known_drug_path=args.known_drug or Path(ct_inputs.get("opentargets_known_drug", "")),
        gene_drug_output=ct_out / "opentargets_gene_drug.parquet",
        drug_indication_output=ct_out / "opentargets_drug_indication.parquet",
    )


if __name__ == "__main__":
    main()
