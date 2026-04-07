#!/usr/bin/env python3
"""Parse TrialPanorama data to standardized gene-drug and drug-indication parquets.

Creates:
1. clinical_trials/trialpanorama_gene_drug.parquet - Gene-drug pairs filtered by pChEMBL threshold
2. clinical_trials/trialpanorama_drug_indication.parquet - Drug-indication with MESH and EFO IDs, numeric phases
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)

# Defaults (can be overridden via config)
DEFAULT_PHASE_STANDARDIZATION = {
    "PHASE1/PHASE2": "PHASE2",
    "PHASE1/PHASE2/PHASE3": "PHASE3",
    "PHASE2/PHASE3": "PHASE3",
    "EARLY_PHASE1": "PHASE1",
}
DEFAULT_VALID_PHASES = ["PHASE1", "PHASE2", "PHASE3", "PHASE4"]
DEFAULT_TRIAL_TYPES = ["INTERVENTIONAL"]
DEFAULT_EXCLUDED_STATUSES = ["", "UNKNOWN"]


def load_trialpanorama_data(input_dir: Path) -> tuple[pd.DataFrame, ...]:
    """Load all TrialPanorama parquet files."""
    log.info(f"Loading TrialPanorama data from {input_dir}")

    drug_moa = pd.read_parquet(input_dir / "drug_moa.parquet")
    log.info(f"Drug MOA: {len(drug_moa):,} records")

    drugs = pd.concat(
        [
            pd.read_parquet(input_dir / "drugs.parquet"),
            pd.read_parquet(input_dir / "drugs_part_2.parquet"),
        ],
        ignore_index=True,
    )
    log.info(f"Drugs: {len(drugs):,} records")

    studies = pd.concat(
        [
            pd.read_parquet(input_dir / "studies.parquet"),
            pd.read_parquet(input_dir / "studies_part_2.parquet"),
        ],
        ignore_index=True,
    )
    log.info(f"Studies: {len(studies):,} records")

    conditions = pd.concat(
        [
            pd.read_parquet(input_dir / "conditions.parquet"),
            pd.read_parquet(input_dir / "conditions_part_2.parquet"),
        ],
        ignore_index=True,
    )
    log.info(f"Conditions: {len(conditions):,} records")

    return drug_moa, drugs, studies, conditions


def create_gene_drug_mapping(
    drug_moa: pd.DataFrame,
    drugs: pd.DataFrame,
    min_pchembl: float,
) -> pd.DataFrame:
    """Create gene-drug mapping filtered by pChEMBL threshold."""
    # Clean drug MOA
    df = drug_moa.dropna(subset=["gene", "drug_name"]).copy()
    df = df[df["organism"] == "Homo sapiens"]
    df = df[~df["gene"].str.contains(r"\|", na=False)]
    df = df[df["act_source"] != "UNKNOWN"]
    log.info(f"After base cleaning: {len(df):,} records")

    # Filter by pChEMBL
    df = df[(df["act_value"].notna()) & (df["act_value"] > min_pchembl)]
    log.info(f"After pChEMBL > {min_pchembl}: {len(df):,} records")

    # Select and rename columns to standard schema
    result = (
        df[["gene", "drug_name", "drug_moa_id", "act_value", "act_source", "action_type"]]
        .rename(columns={"act_value": "pchembl", "act_source": "subsource"})
        .drop_duplicates()
    )

    # Add drugbank_id from drugs table
    if "drugbank_id" in drugs.columns:
        drugbank_map = drugs[["drug_name", "drugbank_id"]].drop_duplicates()
        result = result.merge(drugbank_map, on="drug_name", how="left")

    log.info(f"Gene-drug: {len(result):,} records, {result['gene'].nunique():,} genes")
    return result


def create_drug_trial_mapping(
    drugs: pd.DataFrame,
    studies: pd.DataFrame,
    conditions: pd.DataFrame,
    mesh_to_efo_path: Path,
    *,
    trial_types: list[str] | None = None,
    valid_phases: list[str] | None = None,
    excluded_statuses: list[str] | None = None,
    phase_standardization: dict[str, str] | None = None,
) -> pd.DataFrame:
    """Create drug-trial mapping with indications.

    Args:
        trial_types: Trial types to include (default: ["INTERVENTIONAL"])
        valid_phases: Phases to include after standardization (default: PHASE1-4)
        excluded_statuses: Recruitment statuses to exclude (default: ["", "UNKNOWN"])
        phase_standardization: Mapping to standardize phase names
    """
    trial_types = trial_types or DEFAULT_TRIAL_TYPES
    valid_phases = valid_phases or DEFAULT_VALID_PHASES
    excluded_statuses = excluded_statuses or DEFAULT_EXCLUDED_STATUSES
    phase_standardization = phase_standardization or DEFAULT_PHASE_STANDARDIZATION

    # Standardize phases
    studies = studies.copy()
    studies["phase"] = studies["phase"].replace(phase_standardization)

    # Filter studies
    studies = studies[
        (studies["trial_type"].isin(trial_types))
        & (studies["phase"].isin(valid_phases))
        & (~studies["recruitment_status"].isin(excluded_statuses))
    ]
    log.info(f"Filtered studies: {len(studies):,} (types={trial_types}, phases={valid_phases})")

    # Get study source column
    if "study_source" not in studies.columns:
        studies["study_source"] = studies.get("source", "TRIALPANORAMA")

    # Filter drugs to valid studies
    drugs = drugs[drugs["study_id"].isin(studies["study_id"])]

    # Filter conditions with MESH IDs
    conditions = conditions[conditions["condition_mesh_id"].notna()]

    # Build mapping
    drug_cols = ["study_id", "drug_name", "drug_moa_id"]
    if "drugbank_id" in drugs.columns:
        drug_cols.append("drugbank_id")

    result = (
        drugs[drug_cols]
        .merge(
            studies[["study_id", "phase", "study_source"]],
            on="study_id",
            how="inner",
        )
        .merge(
            conditions[["study_id", "condition_mesh_id", "condition_name"]],
            on="study_id",
            how="left",
        )
    )
    result = result.dropna(subset=["condition_mesh_id"])
    log.info(f"Drug-trial mapping: {len(result):,} records")

    # Map MESH to EFO
    if mesh_to_efo_path.exists():
        efo_map = pd.read_csv(mesh_to_efo_path, sep="\t")
        efo_map["mesh_id"] = efo_map["curie_id"].str.replace("MeSH:", "", regex=False)
        efo_map = efo_map[["mesh_id", "mapped_curie", "mapped_label"]].rename(
            columns={"mapped_curie": "efo_id", "mapped_label": "efo_label"},
        )
        result = result.merge(efo_map, left_on="condition_mesh_id", right_on="mesh_id", how="left")
        result = result.drop(columns=["mesh_id"])
        log.info(
            f"After EFO mapping: {len(result):,} records, {result['efo_id'].notna().sum():,} with EFO",
        )
    else:
        log.warning(f"MESH to EFO mapping not found: {mesh_to_efo_path}")
        result["efo_id"] = None
        result["efo_label"] = None

    # Map string phases to numeric and rename columns to standard schema
    phase_map = {"PHASE1": 1, "PHASE2": 2, "PHASE3": 3, "PHASE4": 4}
    result["phase"] = result["phase"].map(phase_map)
    result = result[result["phase"].notna()].copy()
    result["phase"] = result["phase"].astype("Int64")

    return result.rename(
        columns={
            "condition_mesh_id": "mesh_id",
            "condition_name": "mesh_heading",
            "efo_label": "efo_term",
            "study_source": "subsource",
            "study_id": "study_ids",
        },
    )


def parse_trialpanorama(
    input_dir: Path,
    mesh_to_efo_path: Path,
    gene_drug_output: Path,
    drug_indication_output: Path,
    min_pchembl: float = 7.0,
    *,
    trial_types: list[str] | None = None,
    valid_phases: list[str] | None = None,
    excluded_statuses: list[str] | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse TrialPanorama data."""
    drug_moa, drugs, studies, conditions = load_trialpanorama_data(input_dir)

    gene_drug = create_gene_drug_mapping(drug_moa, drugs, min_pchembl)
    drug_trial = create_drug_trial_mapping(
        drugs,
        studies,
        conditions,
        mesh_to_efo_path,
        trial_types=trial_types,
        valid_phases=valid_phases,
        excluded_statuses=excluded_statuses,
    )

    # Save
    gene_drug_output.parent.mkdir(parents=True, exist_ok=True)
    drug_indication_output.parent.mkdir(parents=True, exist_ok=True)
    gene_drug.to_parquet(gene_drug_output, index=False)
    drug_trial.to_parquet(drug_indication_output, index=False)
    log.info(f"Saved to {gene_drug_output} and {drug_indication_output}")

    return gene_drug, drug_trial


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--input-dir", type=Path, help="TrialPanorama LOCAL_DIR")
    parser.add_argument("--mesh-to-efo", type=Path, help="MESH to EFO mapping TSV")
    parser.add_argument("--output-dir", type=Path, help="Output directory")
    parser.add_argument("--min-pchembl", type=float, help="Minimum pChEMBL threshold")
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
        mappings = cfg.get("mappings", {})
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
        tp_filters = cfg.get("trialpanorama_filters", {})
    else:
        ct_inputs, mappings, thresholds, output_dir, tp_filters = {}, {}, {}, Path(), {}

    ct_out = args.output_dir or output_dir / "clinical_trials"
    parse_trialpanorama(
        input_dir=args.input_dir or Path(ct_inputs.get("trialpanorama_dir", "")),
        mesh_to_efo_path=args.mesh_to_efo or Path(mappings.get("mesh_to_efo", "")),
        gene_drug_output=ct_out / "trialpanorama_gene_drug.parquet",
        drug_indication_output=ct_out / "trialpanorama_drug_indication.parquet",
        min_pchembl=args.min_pchembl
        if args.min_pchembl is not None
        else thresholds.get("chembl_min_pchembl", 7.0),
        trial_types=tp_filters.get("trial_types"),
        valid_phases=tp_filters.get("valid_phases"),
        excluded_statuses=tp_filters.get("excluded_statuses"),
    )


if __name__ == "__main__":
    main()
