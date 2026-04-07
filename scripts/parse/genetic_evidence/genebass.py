#!/usr/bin/env python3
"""Parse Genebass gene-based association results to gene-disease pairs.

Reads:
- genebass_associations: Pre-filtered Genebass parquet (p < 1e-5)
- ukb_phenotype_manifest: UK Biobank phenotype manifest with EFO mappings

Outputs:
- genebass.parquet with columns: gene, efo_id, trait, p_value, source
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def parse_genebass(
    associations_path: Path,
    manifest_path: Path,
    output_path: Path,
    max_pvalue: float = 1e-5,
) -> pd.DataFrame:
    """Parse Genebass to gene-disease associations with EFO mapping."""
    log.info(f"Loading Genebass from {associations_path}")
    df = pd.read_parquet(associations_path)
    log.info(f"Loaded {len(df):,} associations")

    # Filter by p-value
    df = df[df["Pvalue"] <= max_pvalue].copy()
    log.info(f"After p-value <= {max_pvalue}: {len(df):,}")

    # Load phenotype manifest for EFO mapping
    log.info(f"Loading phenotype manifest from {manifest_path}")
    manifest = pd.read_csv(
        manifest_path,
        sep="\t",
        usecols=["phenocode", "coding_description", "trait_efos", "description"],
        dtype=str,
    )

    # Get phenocode -> EFO mapping, exploding comma-separated EFO IDs
    mapping = manifest[manifest["trait_efos"].notna()].copy()
    mapping["efo_id"] = mapping["trait_efos"].str.split(",")
    mapping = mapping.explode("efo_id")
    mapping["efo_id"] = mapping["efo_id"].str.strip()
    mapping = mapping.rename(columns={"description": "manifest_trait"})
    log.info(f"Phenotype manifest has {mapping['phenocode'].nunique():,} phenocodes with EFO")

    # Join on (phenocode, coding_description) to avoid cross-joining coded
    # phenotypes like 20002 ("Non-cancer illness code") across all sub-codings.
    # Quantitative traits have no coding_description (NaN on both sides) — fill
    # with a sentinel so the join still works for those.
    join_cols = ["phenocode", "coding_description"]
    df["phenocode"] = df["phenocode"].astype(str)
    df["coding_description"] = df["coding_description"].fillna("")
    mapping["coding_description"] = mapping["coding_description"].fillna("")
    df = df.merge(
        mapping[["phenocode", "coding_description", "efo_id", "manifest_trait"]],
        on=join_cols,
        how="inner",
    )
    log.info(f"After EFO mapping: {len(df):,} associations")

    # Convert ontology ID format (EFO_xxx -> EFO:xxx, HP_xxx -> HP:xxx, etc.)
    df["efo_id"] = df["efo_id"].str.replace("_", ":", n=1)

    # Use manifest trait if available, else Genebass description
    df["trait"] = df["manifest_trait"].fillna(df["description"])

    # Select and rename columns
    result = df[["gene_symbol", "efo_id", "trait", "Pvalue", "annotation"]].copy()
    result = result.rename(columns={"gene_symbol": "gene", "Pvalue": "p_value"})
    result["source"] = "Genebass"

    # Filter valid genes
    result = result[result["gene"].notna() & (result["gene"] != "")]

    # Keep best p-value per gene-disease pair
    result = result.sort_values("p_value").groupby(["gene", "efo_id"], as_index=False).first()
    log.info(f"Unique gene-disease pairs: {len(result):,}")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved {len(result):,} associations to {output_path}")

    return result


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--input", type=Path, help="Genebass associations parquet")
    parser.add_argument("--manifest", type=Path, help="UK Biobank phenotype manifest TSV")
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
        mappings = cfg.get("mappings", {})
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ge_inputs, mappings, thresholds, output_dir = {}, {}, {}, Path()

    parse_genebass(
        associations_path=args.input or Path(ge_inputs.get("genebass_associations", "")),
        manifest_path=args.manifest or Path(mappings.get("ukb_phenotype_manifest", "")),
        output_path=args.output or output_dir / "genetic_evidence" / "genebass.parquet",
        max_pvalue=args.max_pvalue or thresholds.get("genebass_max_pvalue", 1e-5),
    )


if __name__ == "__main__":
    main()
