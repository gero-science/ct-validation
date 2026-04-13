#!/usr/bin/env python3
"""Parse STITCH protein-chemical interactions to gene-drug mapping."""

import argparse
import logging
from pathlib import Path

import polars as pl
import yaml

log = logging.getLogger(__name__)


def create_ensp_to_hugo_mapping(protein_aliases: pl.DataFrame) -> dict[str, str]:
    """Create ENSP -> HUGO gene symbol mapping."""
    source_priority = {
        "BioMart_HUGO": 1,
        "Ensembl_HGNC_symbol": 2,
        "Ensembl_HGNC": 3,
        "UniProt_GN_Name": 4,
    }

    hugo = protein_aliases.filter(pl.col("source").is_in(list(source_priority.keys())))
    hugo = hugo.with_columns(pl.col("source").replace(source_priority).cast(pl.Int32).alias("priority"))
    hugo = hugo.sort("priority").unique(subset=["#string_protein_id"], keep="first")

    mapping = dict(zip(hugo["#string_protein_id"].to_list(), hugo["alias"].to_list()))
    log.info(f"Mapped {len(mapping):,} STRING protein IDs to HUGO symbols")
    return mapping


def create_chemical_id_mappings(
    chemical_sources: pl.DataFrame,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Create CID -> ChEMBL and CID -> DrugBank mappings."""
    chembl = chemical_sources.filter(pl.col("source_type") == "ChEMBL")
    chembl_mapping = (
        chembl.select(pl.col("chemical").alias("cid"), pl.col("source_id").alias("chembl_id"))
        .unique(subset=["cid"], keep="first")
    )

    drugbank = chemical_sources.filter(pl.col("source_type") == "DrugBank")
    drugbank_mapping = (
        drugbank.select(pl.col("chemical").alias("cid"), pl.col("source_id").alias("drugbank_id"))
        .unique(subset=["cid"], keep="first")
    )

    log.info(f"ChEMBL mappings: {len(chembl_mapping):,}, DrugBank: {len(drugbank_mapping):,}")
    return chembl_mapping, drugbank_mapping


def parse_stitch(
    actions_path: Path,
    chemical_sources_path: Path,
    protein_aliases_path: Path,
    output_path: Path,
    min_score: int = 700,
) -> pl.DataFrame:
    """Parse STITCH data to gene-drug mapping."""
    log.info("Loading STITCH data...")
    actions = pl.read_csv(actions_path, separator="\t")
    chemical_sources = (
        pl.scan_csv(
            chemical_sources_path,
            separator="\t",
            has_header=False,
            new_columns=["chemical", "alias", "source_type", "source_id"],
            skip_rows=1,
            comment_prefix="#",
        )
        .filter(pl.col("source_type").is_in(["ChEMBL", "DrugBank"]))
        .collect()
    )
    protein_aliases = pl.read_csv(protein_aliases_path, separator="\t")

    ensp_to_hugo = create_ensp_to_hugo_mapping(protein_aliases)
    chembl_mapping, drugbank_mapping = create_chemical_id_mappings(chemical_sources)

    # Filter actions
    log.info(f"Filtering actions (score >= {min_score})...")
    clean = actions.filter(
        pl.col("action").is_in(["activation", "inhibition"])
        & (pl.col("score") >= min_score)
        & pl.col("item_id_a").str.starts_with("CID")
        & pl.col("item_id_b").str.starts_with("9606.")
        & (pl.col("a_is_acting") == "t")
    )

    # Map proteins to genes
    clean = clean.with_columns(
        pl.col("item_id_b").replace_strict(ensp_to_hugo, default=None).alias("gene")
    ).filter(pl.col("gene").is_not_null())
    log.info(f"After gene mapping: {len(clean):,} interactions")

    # Add chemical ID mappings
    clean = clean.join(chembl_mapping, left_on="item_id_a", right_on="cid", how="left")
    clean = clean.join(drugbank_mapping, left_on="item_id_a", right_on="cid", how="left")

    # Select and rename columns
    result = clean.select(
        "gene",
        "chembl_id",
        "drugbank_id",
        pl.col("item_id_a").alias("stitch_id"),
        pl.col("action").alias("action_type"),
        pl.col("mode").alias("interaction_mode"),
        pl.col("score").alias("stitch_score"),
        pl.lit("STITCH").alias("subsource"),
        pl.lit(None).cast(pl.Utf8).alias("drug_name"),
    )

    log.info(f"Result: {len(result):,} interactions, {result['gene'].n_unique():,} genes")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.write_parquet(output_path)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--actions", type=Path, help="STITCH actions.tsv.gz")
    parser.add_argument("--chemical-sources", type=Path, help="STITCH chemical_sources.tsv.gz")
    parser.add_argument("--protein-aliases", type=Path, help="STRING protein_aliases.txt.gz")
    parser.add_argument("--output", type=Path, help="Output parquet path")
    parser.add_argument("--min-score", type=int, help="Minimum STITCH score")
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
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ct_inputs, thresholds, output_dir = {}, {}, Path()

    actions = args.actions or Path(ct_inputs.get("stitch_actions", ""))
    chemical_sources = args.chemical_sources or Path(ct_inputs.get("stitch_chemical_sources", ""))
    protein_aliases = args.protein_aliases or Path(ct_inputs.get("stitch_protein_aliases", ""))
    output = args.output or output_dir / "clinical_trials" / "stitch_gene_drug.parquet"
    min_score = (
        args.min_score if args.min_score is not None else thresholds.get("stitch_min_score", 700)
    )

    parse_stitch(actions, chemical_sources, protein_aliases, output, min_score)


if __name__ == "__main__":
    main()
