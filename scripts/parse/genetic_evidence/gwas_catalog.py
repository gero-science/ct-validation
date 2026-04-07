#!/usr/bin/env python3
"""Parse GWAS Catalog to gene-trait associations."""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def split_genes(gene_str) -> list[str]:
    """Split gene string by common delimiters."""
    if pd.isna(gene_str):
        return []
    genes = str(gene_str).replace(" - ", ",").replace(" x ", ",").split(",")
    return [g.strip() for g in genes if g.strip()]


def extract_ontology_id(uri) -> str | None:
    """Extract ontology ID from URI (e.g., '.../EFO_0000712' -> 'EFO:0000712')."""
    if pd.isna(uri):
        return None
    parts = str(uri).split("/")[-1]
    if "_" in parts:
        return parts.replace("_", ":", 1)
    return parts


def parse_gwas_catalog(
    raw_path: Path,
    output_path: Path,
    max_pvalue: float = 1e-8,
) -> pd.DataFrame:
    """Parse GWAS Catalog to gene-trait associations."""
    log.info("Loading GWAS Catalog...")
    df = pd.read_csv(raw_path, sep="\t", low_memory=False)
    log.info(f"Loaded {len(df):,} associations")

    # Filter by p-value
    df = df[df["P-VALUE"] < max_pvalue].copy()
    log.info(f"After p-value < {max_pvalue}: {len(df):,}")

    # Filter for mapped genes and traits
    df = df[df["MAPPED_GENE"].notna() & df["MAPPED_TRAIT"].notna()].copy()
    log.info(f"After mapped gene/trait filter: {len(df):,}")

    # Filter out multi-EFO traits
    df["n_uris"] = df["MAPPED_TRAIT_URI"].apply(
        lambda x: len(str(x).split(",")) if pd.notna(x) else 0,
    )
    df = df[df["n_uris"] == 1].copy()
    log.info(f"After single-EFO filter: {len(df):,}")

    # Split and explode genes
    df["GENES"] = df["MAPPED_GENE"].apply(split_genes)
    df = df.explode("GENES").copy()
    log.info(f"After gene explosion: {len(df):,}")

    # Create gene-trait pairs with p-value preserved
    pairs = df[["GENES", "MAPPED_TRAIT", "MAPPED_TRAIT_URI", "P-VALUE"]].copy()
    pairs.columns = ["gene", "trait", "trait_uri", "p_value"]

    # Clean up
    pairs = pairs[pairs["gene"].notna() & (pairs["gene"] != "") & pairs["trait"].notna()].copy()
    pairs["gene"] = pairs["gene"].str.upper()

    # Extract ontology ID
    pairs["efo_id"] = pairs["trait_uri"].apply(extract_ontology_id)

    # Deduplicate (keep lowest p-value)
    pairs = pairs.sort_values("p_value").groupby(["gene", "efo_id"], as_index=False).first()
    log.info(f"Unique gene-disease pairs: {len(pairs):,}")

    # Final output - keep p_value for downstream filtering
    pairs["source"] = "GWAS"
    result = pairs[["gene", "efo_id", "trait", "p_value", "source"]].copy()

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved {len(result):,} associations to {output_path}")

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--input", type=Path, help="GWAS Catalog TSV (gzipped)")
    parser.add_argument("--output", type=Path, help="Output parquet path")
    parser.add_argument("--max-pvalue", type=float, help="Maximum p-value threshold")
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
        ge_inputs = cfg.get("genetic_evidence", {})
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ge_inputs, thresholds, output_dir = {}, {}, Path()

    parse_gwas_catalog(
        raw_path=args.input or Path(ge_inputs.get("gwas_catalog", "")),
        output_path=args.output or output_dir / "genetic_evidence" / "gwas_catalog.parquet",
        max_pvalue=args.max_pvalue
        if args.max_pvalue is not None
        else thresholds.get("gwas_max_pvalue", 1e-8),
    )


if __name__ == "__main__":
    main()
