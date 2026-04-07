#!/usr/bin/env python3
"""Parse STITCH protein-chemical interactions to gene-drug mapping."""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def create_ensp_to_hugo_mapping(protein_aliases: pd.DataFrame) -> dict[str, str]:
    """Create ENSP -> HUGO gene symbol mapping."""
    hugo_sources = protein_aliases[
        protein_aliases["source"].isin(
            ["BioMart_HUGO", "Ensembl_HGNC_symbol", "Ensembl_HGNC", "UniProt_GN_Name"],
        )
    ].copy()

    source_priority = {
        "BioMart_HUGO": 1,
        "Ensembl_HGNC_symbol": 2,
        "Ensembl_HGNC": 3,
        "UniProt_GN_Name": 4,
    }
    hugo_sources["priority"] = hugo_sources["source"].map(source_priority)
    hugo_sources = hugo_sources.sort_values("priority").drop_duplicates(
        subset=["#string_protein_id"],
        keep="first",
    )

    mapping = dict(zip(hugo_sources["#string_protein_id"], hugo_sources["alias"]))
    log.info(f"Mapped {len(mapping):,} STRING protein IDs to HUGO symbols")
    return mapping


def create_chemical_id_mappings(
    chemical_sources: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create CID -> ChEMBL and CID -> DrugBank mappings."""
    # ChEMBL mapping
    chembl_sources = chemical_sources[chemical_sources["source_type"] == "ChEMBL"][
        ["chemical", "alias", "source_id"]
    ]
    chembl_from_chemical = chembl_sources[["chemical", "source_id"]].rename(
        columns={"chemical": "cid", "source_id": "chembl_id"},
    )
    chembl_from_alias = chembl_sources[["alias", "source_id"]].rename(
        columns={"alias": "cid", "source_id": "chembl_id"},
    )
    chembl_mapping = pd.concat([chembl_from_chemical, chembl_from_alias], ignore_index=True)
    chembl_mapping["chembl_id"] = chembl_mapping["chembl_id"].astype(str)
    chembl_mapping = chembl_mapping.drop_duplicates(subset=["cid"], keep="first")

    # DrugBank mapping
    drugbank_sources = chemical_sources[chemical_sources["source_type"] == "DrugBank"][
        ["chemical", "alias", "source_id"]
    ]
    drugbank_from_chemical = drugbank_sources[["chemical", "source_id"]].rename(
        columns={"chemical": "cid", "source_id": "drugbank_id"},
    )
    drugbank_from_alias = drugbank_sources[["alias", "source_id"]].rename(
        columns={"alias": "cid", "source_id": "drugbank_id"},
    )
    drugbank_mapping = pd.concat([drugbank_from_chemical, drugbank_from_alias], ignore_index=True)
    drugbank_mapping = drugbank_mapping.drop_duplicates(subset=["cid"], keep="first")

    log.info(f"ChEMBL mappings: {len(chembl_mapping):,}, DrugBank: {len(drugbank_mapping):,}")
    return chembl_mapping, drugbank_mapping


def parse_stitch(
    actions_path: Path,
    chemical_sources_path: Path,
    protein_aliases_path: Path,
    output_path: Path,
    min_score: int = 700,
) -> pd.DataFrame:
    """Parse STITCH data to gene-drug mapping."""
    log.info("Loading STITCH data...")
    actions = pd.read_csv(actions_path, sep="\t")
    chemical_sources = pd.read_csv(
        chemical_sources_path,
        sep="\t",
        comment="#",
        names=["chemical", "alias", "source_type", "source_id"],
        skiprows=1,
        low_memory=False,
    )
    protein_aliases = pd.read_csv(protein_aliases_path, sep="\t")

    # Create mappings
    ensp_to_hugo = create_ensp_to_hugo_mapping(protein_aliases)
    chembl_mapping, drugbank_mapping = create_chemical_id_mappings(chemical_sources)

    # Filter actions
    log.info(f"Filtering actions (score >= {min_score})...")
    clean_actions = actions[
        actions["action"].isin(["activation", "inhibition"])
        & (actions["score"] >= min_score)
        & actions["item_id_a"].str.startswith("CID")
        & actions["item_id_b"].str.startswith("9606.")
        & (actions["a_is_acting"] == "t")
    ].copy()

    # Map proteins to genes
    clean_actions["gene"] = clean_actions["item_id_b"].map(ensp_to_hugo)
    clean_actions = clean_actions[clean_actions["gene"].notna()].copy()
    log.info(f"After gene mapping: {len(clean_actions):,} interactions")

    # Add chemical ID mappings
    clean_actions = clean_actions.merge(
        chembl_mapping,
        left_on="item_id_a",
        right_on="cid",
        how="left",
    )
    clean_actions = clean_actions.drop(columns=["cid"], errors="ignore")
    clean_actions = clean_actions.merge(
        drugbank_mapping,
        left_on="item_id_a",
        right_on="cid",
        how="left",
    )
    clean_actions = clean_actions.drop(columns=["cid"], errors="ignore")

    # Select and rename columns
    result = clean_actions[
        ["gene", "chembl_id", "drugbank_id", "item_id_a", "action", "mode", "score"]
    ].rename(
        columns={
            "item_id_a": "stitch_id",
            "action": "action_type",
            "mode": "interaction_mode",
            "score": "stitch_score",
        },
    )

    result["subsource"] = "STITCH"
    result["drug_name"] = None

    log.info(f"Result: {len(result):,} interactions, {result['gene'].nunique():,} genes")

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")

    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--config",
        type=Path,
        help="Config YAML (uses configs/parsing.yaml by default)",
    )
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

    # Load config
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

    # CLI args override config
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
