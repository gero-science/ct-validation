#!/usr/bin/env python3
"""Aggregate genetic evidence sources into unified gene-disease pairs."""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def load_and_normalize(
    path: Path,
    source: str,
    efo_col: str = "efo_id",
    pvalue_col: str | None = None,
    score_col: str | None = None,
) -> pd.DataFrame:
    """Load parquet and normalize to standard columns."""
    df = pd.read_parquet(path)
    if efo_col != "efo_id":
        df = df.rename(columns={efo_col: "efo_id"})
    if "trait" not in df.columns:
        df["trait"] = None
    df["source"] = source
    # Add score columns
    df["gwas_pvalue"] = df[pvalue_col] if pvalue_col and pvalue_col in df.columns else None
    df["ot_score"] = df[score_col] if score_col and score_col in df.columns else None
    return df[["gene", "efo_id", "trait", "source", "gwas_pvalue", "ot_score"]]


def aggregate_genetic_evidence(
    gwas_path: Path,
    clinvar_path: Path,
    opentargets_path: Path,
    genebass_path: Path,
    output_path: Path,
) -> pd.DataFrame:
    """Aggregate genetic evidence sources."""
    # Load sources with score columns where available
    sources = []
    source_configs = [
        (gwas_path, "GWAS Catalog", "efo_id", "p_value", None),
        (clinvar_path, "ClinVar", "efo_id", None, None),
        (opentargets_path, "OpenTargets", "efo_id", None, "score"),
        (genebass_path, "Genebass", "efo_id", "p_value", None),
    ]
    for path, name, efo_col, pval_col, score_col in source_configs:
        if path.exists():
            df = load_and_normalize(path, name, efo_col, pval_col, score_col)
            sources.append(df)
            log.info(f"Loaded {name}: {len(df):,} pairs")
        else:
            log.warning(f"Missing: {path}")

    # Concatenate and aggregate
    combined = pd.concat(sources, ignore_index=True)
    log.info(f"Combined: {len(combined):,} total records")

    result = combined.groupby(["gene", "efo_id"], as_index=False).agg(
        {
            "trait": "first",
            "source": lambda x: "&".join(sorted(set(x))),
            "gwas_pvalue": "min",  # Best (lowest) p-value
            "ot_score": "max",  # Best (highest) OT score
        },
    )
    log.info(f"Aggregated: {len(result):,} unique gene-disease pairs")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--gwas", type=Path, help="GWAS Catalog parquet")
    parser.add_argument("--clinvar", type=Path, help="ClinVar parquet")
    parser.add_argument("--opentargets", type=Path, help="OpenTargets parquet")
    parser.add_argument("--genebass", type=Path, help="Genebass parquet")
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
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        output_dir = Path()

    ge_dir = output_dir / "genetic_evidence"
    aggregate_genetic_evidence(
        gwas_path=args.gwas or ge_dir / "gwas_catalog.parquet",
        clinvar_path=args.clinvar or ge_dir / "clinvar.parquet",
        opentargets_path=args.opentargets or ge_dir / "opentargets.parquet",
        genebass_path=args.genebass or ge_dir / "genebass.parquet",
        output_path=args.output or ge_dir / "aggregated" / "genetic_evidence.parquet",
    )


if __name__ == "__main__":
    main()
