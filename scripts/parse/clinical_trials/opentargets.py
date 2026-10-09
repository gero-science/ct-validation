#!/usr/bin/env python3
"""Parse Open Targets clinical data to standardized gene-drug and drug-indication parquets.

Reads, up to 25.12:
- opentargets_known_drug: Open Targets known_drug parquet directory

or, from 26.06, the datasets that replace it:
- opentargets_clinical_report: drug-indication reports
- opentargets_drug_mechanism_of_action: drug-target links
- opentargets_drug_molecule: drug names and metadata
- opentargets_target: target symbols

26.03, the first release without `known_drug`, is not supported: it aggregated its reports
under another rule.

Outputs:
- clinical_trials/opentargets_gene_drug.parquet
- clinical_trials/opentargets_drug_indication.parquet
"""

import argparse
import logging
from pathlib import Path

import duckdb
import pandas as pd
import trial_status
import yaml
from ct_validation.efo import EfoOntology, efo_for_remap, remap_efo_ids

log = logging.getLogger(__name__)

# Clinical stages from 26.06 on the internal 1-4 scale. Combined phases take their highest
# component, as TrialPanorama's do, and a market withdrawal follows an approval. IND,
# PRECLINICAL and UNKNOWN establish no clinical phase and are dropped.
CLINICAL_STAGE_PHASE = {
    "EARLY_PHASE_1": 1,
    "PHASE_1": 1,
    "PHASE_1_2": 2,
    "PHASE_2": 2,
    "PHASE_2_3": 3,
    "PHASE_3": 3,
    "PREAPPROVAL": 4,
    "PHASE_4": 4,
    "APPROVAL": 4,
    "WITHDRAWAL": 4,
}

# The quality controls that keep a report out of Open Targets' own `clinical_indication`.
# UNVALIDATED_INDICATION is a warning there rather than a filter, and NO_DISEASE comes with
# an unmapped disease, which is dropped anyway.
EXCLUDED_QUALITY_CONTROLS = {"INDIRECT_PRIMARY_PURPOSE", "PHASE_IV_NOT_APPROVED"}


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
        ["drugId", "prefName", "drugType", "diseaseId", "label", "phase", "status", "urls"]
    ].copy()

    # A third of the ids are MONDO, HP, Orphanet, MP, GO or DOID. EFO is a merged ontology
    # that imports those under their own prefixes, so keep the id as recorded rather than
    # rewriting it to an EFO ancestor.
    result["efo_id"] = result["diseaseId"].str.replace("_", ":", n=1)

    # Filter to rows with valid EFO and phase
    result["phase"] = pd.to_numeric(result["phase"], errors="coerce")
    result = result[(result["phase"].notna()) & (result["phase"] >= 1) & (result["phase"] <= 4)]
    result = result[result["phase"] == result["phase"].astype(int)]  # drop non-integer (0.5)
    result["phase"] = result["phase"].astype("Int64")

    # Extract subsources from urls
    result["subsource"] = result["urls"].apply(_extract_subsources)

    result = trial_status.apply(
        result,
        trial_status.OPENTARGETS_CONCLUDED,
        trial_status.OPENTARGETS_DROPPED,
    )

    result = result.rename(
        columns={
            "drugId": "chembl_id",
            "prefName": "drug_name",
            "drugType": "molecule_type",
            "label": "efo_term",
        },
    )

    keep = [
        "efo_id",
        "phase",
        "status",
        "is_concluded",
        "drug_name",
        "chembl_id",
        "efo_term",
        "subsource",
        "molecule_type",
    ]
    result = result[keep].drop_duplicates()

    log.info(
        f"Drug-indication: {len(result):,} rows, {result['chembl_id'].nunique():,} drugs, "
        f"{result['efo_id'].nunique():,} diseases",
    )
    return result


def _molecules_sql(drug_molecule_path: Path) -> str:
    """SQL for drug names and metadata from Open Targets drug_molecule, keyed by ChEMBL id."""
    return f"""
        SELECT id AS chembl_id,
               name AS drug_name,
               drugType AS molecule_type,
               NULLIF(array_to_string([t.label FOR t IN tradeNames], '|'), '') AS trade_names,
               NULLIF(array_to_string([s.label FOR s IN synonyms], '|'), '') AS synonyms
        FROM read_parquet('{drug_molecule_path / "*.parquet"}')
    """


def parse_mechanisms(
    drug_mechanism_of_action_path: Path,
    drug_molecule_path: Path,
    target_path: Path,
) -> pd.DataFrame:
    """Extract gene-drug pairs from Open Targets drug_mechanism_of_action (26.06 onwards)."""
    with duckdb.connect() as con:
        result = con.sql(
            # sql
            f"""
            WITH links AS (
                SELECT chembl_id, target_id, actionType AS action_type
                FROM read_parquet('{drug_mechanism_of_action_path / "*.parquet"}'),
                     unnest(chemblIds) AS c(chembl_id),
                     unnest(targets) AS t(target_id)
            )
            SELECT DISTINCT t.approvedSymbol AS gene, m.drug_name, l.chembl_id,
                   m.molecule_type, l.action_type, m.trade_names, m.synonyms
            FROM links l
            JOIN read_parquet('{target_path / "*.parquet"}') t ON t.id = l.target_id
            LEFT JOIN ({_molecules_sql(drug_molecule_path)}) m USING (chembl_id)
            WHERE t.approvedSymbol IS NOT NULL AND l.chembl_id IS NOT NULL
            """,
        ).df()

    log.info(f"Gene-drug: {len(result):,} rows, {result['gene'].nunique():,} genes")
    return result


def parse_clinical_reports(clinical_report_path: Path, drug_molecule_path: Path) -> pd.DataFrame:
    """Extract drug-indication records from Open Targets clinical_report (26.06 onwards).

    One record per report, drug and disease, under the quality controls Open Targets
    applies when it aggregates them into `clinical_indication`. Registry trials are kept
    only when interventional; labels, agencies and curated resources carry no study type.
    `type` is not read: it holds the report's origin up to 26.06 and INDICATION or SAFETY
    from 26.09, and safety reports name no disease in either, so the unnest drops them.
    """
    phase_cases = " ".join(f"WHEN '{s}' THEN {p}" for s, p in CLINICAL_STAGE_PHASE.items())
    excluded = ", ".join(f"'{qc}'" for qc in sorted(EXCLUDED_QUALITY_CONTROLS))
    with duckdb.connect() as con:
        result = con.sql(
            # sql
            f"""
            WITH records AS (
                SELECT drug.drugId AS chembl_id,
                       -- First underscore only, matching the GE side these pairs join against.
                       regexp_replace(disease.diseaseId, '_', ':') AS efo_id,
                       CASE clinicalStage {phase_cases} END AS phase,
                       trialOverallStatus AS status,
                       source AS subsource
                FROM read_parquet('{clinical_report_path / "*.parquet"}'),
                     unnest(drugs) AS u(drug),
                     unnest(diseases) AS v(disease)
                WHERE NOT list_has_any(COALESCE(qualityControls, []), [{excluded}])
                  AND COALESCE(trialStudyType, 'INTERVENTIONAL') = 'INTERVENTIONAL'
            )
            SELECT DISTINCT r.efo_id, r.phase, r.status, m.drug_name, r.chembl_id,
                   r.subsource, m.molecule_type
            FROM records r
            LEFT JOIN ({_molecules_sql(drug_molecule_path)}) m USING (chembl_id)
            WHERE r.chembl_id IS NOT NULL AND r.efo_id IS NOT NULL AND r.phase IS NOT NULL
            """,
        ).df()

    result = trial_status.apply(
        result,
        trial_status.OPENTARGETS_REPORT_CONCLUDED,
        trial_status.OPENTARGETS_REPORT_DROPPED,
    )

    log.info(
        f"Drug-indication: {len(result):,} rows, {result['chembl_id'].nunique():,} drugs, "
        f"{result['efo_id'].nunique():,} diseases",
    )
    return result


def parse_opentargets(
    gene_drug_output: Path,
    drug_indication_output: Path,
    known_drug_path: Path | None = None,
    clinical_report_path: Path | None = None,
    drug_mechanism_of_action_path: Path | None = None,
    drug_molecule_path: Path | None = None,
    target_path: Path | None = None,
    efo: EfoOntology | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse Open Targets to standardized outputs.

    Takes either `known_drug_path` (up to 25.12) or the four datasets that replace it from
    26.06: `clinical_report`, `drug_mechanism_of_action`, `drug_molecule` and `target`.
    With `efo`, obsolete disease ids are rewritten to their successor.
    """
    if known_drug_path is not None:
        log.info(f"Loading Open Targets known_drug from {known_drug_path}")
        df = pd.read_parquet(known_drug_path)
        gene_drug = parse_gene_drug(df)
        drug_indication = parse_drug_indication(df)
        if efo is not None:
            drug_indication = remap_efo_ids(drug_indication, efo, label_column="efo_term")
    elif (
        clinical_report_path is not None
        and drug_mechanism_of_action_path is not None
        and drug_molecule_path is not None
        and target_path is not None
    ):
        log.info(f"Loading Open Targets clinical reports from {clinical_report_path}")
        gene_drug = parse_mechanisms(drug_mechanism_of_action_path, drug_molecule_path, target_path)
        drug_indication = parse_clinical_reports(clinical_report_path, drug_molecule_path)
        if efo is not None:
            drug_indication = remap_efo_ids(drug_indication, efo)
    else:
        msg = (
            "Pass either known_drug (up to 25.12) or all of clinical_report, "
            "drug_mechanism_of_action, drug_molecule and target (26.06 onwards)"
        )
        raise ValueError(msg)

    gene_drug_output.parent.mkdir(parents=True, exist_ok=True)
    drug_indication_output.parent.mkdir(parents=True, exist_ok=True)
    gene_drug.to_parquet(gene_drug_output, index=False)
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
    parser.add_argument("--clinical-report", type=Path, help="clinical_report dir (26.06+)")
    parser.add_argument(
        "--drug-mechanism-of-action",
        type=Path,
        help="drug_mechanism_of_action dir (26.06+)",
    )
    parser.add_argument("--drug-molecule", type=Path, help="drug_molecule dir (26.06+)")
    parser.add_argument("--target", type=Path, help="target dir (26.06+)")
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
        mappings = cfg.get("mappings", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ct_inputs, mappings, output_dir = {}, {}, Path()

    def dataset(arg: Path | None, key: str) -> Path | None:
        return arg or (Path(ct_inputs[key]) if ct_inputs.get(key) else None)

    ct_out = args.output_dir or output_dir / "clinical_trials"
    parse_opentargets(
        gene_drug_output=ct_out / "opentargets_gene_drug.parquet",
        drug_indication_output=ct_out / "opentargets_drug_indication.parquet",
        known_drug_path=dataset(args.known_drug, "opentargets_known_drug"),
        clinical_report_path=dataset(args.clinical_report, "opentargets_clinical_report"),
        drug_mechanism_of_action_path=dataset(
            args.drug_mechanism_of_action,
            "opentargets_drug_mechanism_of_action",
        ),
        drug_molecule_path=dataset(args.drug_molecule, "opentargets_drug_molecule"),
        target_path=dataset(args.target, "opentargets_target"),
        efo=efo_for_remap(mappings),
    )


if __name__ == "__main__":
    main()
