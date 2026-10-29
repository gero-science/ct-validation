"""Tests for clinical_trials/opentargets.py's efo_id derivation.

EFO is a merged ontology: MONDO, HP and Orphanet terms are imported under their own CURIE
prefixes. `parse_drug_indication` keeps each disease id as recorded, converting `_` to `:`
and nothing else.
"""

from pathlib import Path

import opentargets
import pandas as pd
import pytest


def _known_drug_row(disease_id, **overrides) -> pd.DataFrame:
    """A single known_drug row with the columns parse_drug_indication reads."""
    row = {
        "drugId": "CHEMBL1",
        "prefName": "Drug1",
        "drugType": "Small molecule",
        "diseaseId": disease_id,
        "label": "some disease",
        "phase": 2,
        "status": "Completed",
        "urls": None,
    }
    row.update(overrides)
    return pd.DataFrame([row])


@pytest.mark.parametrize(
    ("disease_id", "expected"),
    [
        ("EFO_0000616", "EFO:0000616"),
        ("MONDO_0004976", "MONDO:0004976"),
        ("HP_0003326", "HP:0003326"),
        ("Orphanet_100992", "Orphanet:100992"),
    ],
)
def test_disease_id_is_kept_as_recorded(disease_id, expected):
    """Whatever ontology the term comes from, it survives as itself."""
    result = opentargets.parse_drug_indication(_known_drug_row(disease_id))

    assert result["efo_id"].tolist() == [expected]


def test_only_the_first_underscore_becomes_a_colon():
    """A CURIE separates prefix from local id once; later underscores belong to the id."""
    result = opentargets.parse_drug_indication(_known_drug_row("MONDO_0004_976"))

    assert result["efo_id"].tolist() == ["MONDO:0004_976"]


def test_missing_disease_id_yields_no_term():
    """A row with no disease id carries none downstream, where the aggregate drops it."""
    result = opentargets.parse_drug_indication(_known_drug_row(None))

    assert result["efo_id"].isna().all()


# --- From 26.06: clinical_report and drug_mechanism_of_action ---


def _report(**overrides) -> dict:
    """An interventional ClinicalTrials.gov report with one drug and one disease."""
    row = {
        "id": "nct00000001",
        "type": "INDICATION",
        "source": "ClinicalTrials.gov",
        "clinicalStage": "PHASE_2",
        "phaseFromSource": "PHASE2",
        "trialOverallStatus": "COMPLETED",
        "trialStudyType": "INTERVENTIONAL",
        "qualityControls": None,
        "drugs": [{"drugFromSource": "drug1", "drugId": "CHEMBL1"}],
        "diseases": [{"diseaseFromSource": "asthma", "diseaseId": "MONDO_0004979"}],
    }
    row.update(overrides)
    return row


def _label() -> dict:
    """A drug-label report: approval stage, no trial status or study type."""
    return _report(
        id="label1",
        source="DailyMed",
        clinicalStage="APPROVAL",
        phaseFromSource="approval",
        trialOverallStatus=None,
        trialStudyType=None,
    )


def _mechanism(**overrides) -> dict:
    """A drug_mechanism_of_action row linking CHEMBL1 to ENSG01."""
    row = {
        "chemblIds": ["CHEMBL1"],
        "targets": ["ENSG01"],
        "actionType": "INHIBITOR",
        "references": [{"source": "PubMed", "ids": ["1"], "urls": ["https://pubmed/1"]}],
    }
    row.update(overrides)
    return row


def _release(
    tmp_path: Path,
    reports: list[dict],
    mechanisms: list[dict] | None = None,
) -> dict[str, Path]:
    """Write a minimal 26.06 release, one parquet per dataset, keyed by argument name."""
    datasets = {
        "clinical_report": pd.DataFrame(reports),
        "drug_mechanism_of_action": pd.DataFrame(
            mechanisms or [_mechanism()],
        ),
        "drug_molecule": pd.DataFrame(
            [
                {
                    "id": "CHEMBL1",
                    "name": "DRUG1",
                    "drugType": "Small molecule",
                    "synonyms": [{"label": "DRUG-1", "source": "ChEMBL"}],
                    "tradeNames": [{"label": "Brand1", "source": "ChEMBL"}],
                },
            ],
        ),
        "target": pd.DataFrame(
            {
                "id": ["ENSG01", "ENSG02"],
                "approvedSymbol": ["GENEA", "GENEB"],
                "targetClass": [
                    [
                        {"id": 1, "label": "Enzyme", "level": "l1"},
                        {"id": 2, "label": "Kinase", "level": "l2"},
                    ],
                    [],
                ],
            },
        ),
    }
    paths = {}
    for name, df in datasets.items():
        (tmp_path / name).mkdir()
        df.to_parquet(tmp_path / name / "00000000.parquet", index=False)
        paths[f"{name}_path"] = tmp_path / name
    return paths


def _parse_reports(release: dict[str, Path]) -> pd.DataFrame:
    return opentargets.parse_clinical_reports(
        release["clinical_report_path"],
        release["drug_molecule_path"],
    )


@pytest.mark.parametrize(
    ("stage", "phase"),
    [
        ("EARLY_PHASE_1", 1),
        ("PHASE_1_2", 2),
        ("PHASE_2_3", 3),
        ("PREAPPROVAL", 4),
        ("PHASE_4", 4),
        ("WITHDRAWAL", 4),
    ],
)
def test_stages_round_up_to_the_internal_scale(tmp_path, stage, phase):
    """Combined phases take their highest component; a market withdrawal was an approval."""
    release = _release(tmp_path, [_report(clinicalStage=stage)])

    assert _parse_reports(release)["phase"].tolist() == [phase]


@pytest.mark.parametrize("stage", ["IND", "PRECLINICAL", "UNKNOWN"])
def test_stages_without_a_clinical_phase_are_dropped(tmp_path, stage):
    """IND, PRECLINICAL and UNKNOWN show no clinical phase was reached."""
    release = _release(tmp_path, [_report(clinicalStage=stage)])

    assert _parse_reports(release).empty


@pytest.mark.parametrize(
    ("flags", "kept"),
    [
        (["INDIRECT_PRIMARY_PURPOSE"], False),
        (["PHASE_IV_NOT_APPROVED"], False),
        (["UNVALIDATED_INDICATION"], True),
    ],
)
def test_quality_controls_match_clinical_indication(tmp_path, flags, kept):
    """Only the flags that keep a report out of clinical_indication drop it here."""
    release = _release(tmp_path, [_report(qualityControls=flags)])

    assert len(_parse_reports(release)) == int(kept)


@pytest.mark.parametrize("study_type", ["OBSERVATIONAL", "EXPANDED_ACCESS"])
def test_non_interventional_trials_are_dropped(tmp_path, study_type):
    """Registry trials count only when interventional, as TrialPanorama's do."""
    release = _release(tmp_path, [_report(trialStudyType=study_type), _label()])

    assert _parse_reports(release)["subsource"].tolist() == ["DailyMed"]


def test_report_type_is_not_read(tmp_path):
    """Up to 26.06 `type` holds the origin rather than INDICATION, and the report still counts."""
    release = _release(tmp_path, [_report(type="CLINICAL_TRIAL", source="AACT")])

    assert _parse_reports(release)["subsource"].tolist() == ["AACT"]


def test_safety_reports_name_no_disease_and_are_dropped(tmp_path):
    """A safety report carries no disease, in either layout, so it yields no record."""
    safety = _report(id="safety1", type="SAFETY", clinicalStage="WITHDRAWAL", diseases=None)
    release = _release(tmp_path, [_report(), safety])

    assert _parse_reports(release)["phase"].tolist() == [2]


def test_withdrawn_trials_are_dropped(tmp_path):
    """A trial withdrawn before enrolment does not show the phase was reached."""
    release = _release(tmp_path, [_report(trialOverallStatus="WITHDRAWN")])

    assert _parse_reports(release).empty


def test_trial_status_sets_is_concluded(tmp_path):
    """A running trial leaves its phase open; a label is a concluded approval."""
    release = _release(tmp_path, [_report(trialOverallStatus="RECRUITING"), _label()])

    result = _parse_reports(release).sort_values("phase")

    assert result[["phase", "is_concluded"]].values.tolist() == [[2, False], [4, True]]


@pytest.mark.parametrize(
    ("stage", "phase_from_source", "concluded"),
    [
        ("PHASE_2", "phase 2", False),
        ("PHASE_2", "discontinued in phase 2", True),
        ("WITHDRAWAL", "withdrawn", True),
    ],
)
def test_status_less_reports_conclude_only_when_settled(
    tmp_path,
    stage,
    phase_from_source,
    concluded,
):
    """A curated 'phase 2' may still be running; a discontinuation or withdrawal is not."""
    report = _report(
        source="TTD",
        clinicalStage=stage,
        phaseFromSource=phase_from_source,
        trialOverallStatus=None,
        trialStudyType=None,
    )

    result = _parse_reports(_release(tmp_path, [report]))

    assert result["is_concluded"].tolist() == [concluded]


def test_one_record_per_drug_and_disease(tmp_path):
    """A report naming two drugs and two diseases yields four records; unmapped ones none."""
    report = _report(
        drugs=[
            {"drugFromSource": "drug1", "drugId": "CHEMBL1"},
            {"drugFromSource": "drug2", "drugId": "CHEMBL2"},
        ],
        diseases=[
            {"diseaseFromSource": "asthma", "diseaseId": "MONDO_0004979"},
            {"diseaseFromSource": "copd", "diseaseId": "EFO_0000341"},
            {"diseaseFromSource": "wheeze", "diseaseId": None},
        ],
    )

    result = _parse_reports(_release(tmp_path, [report]))

    assert sorted(zip(result["chembl_id"], result["efo_id"], strict=True)) == [
        ("CHEMBL1", "EFO:0000341"),
        ("CHEMBL1", "MONDO:0004979"),
        ("CHEMBL2", "EFO:0000341"),
        ("CHEMBL2", "MONDO:0004979"),
    ]


def test_mechanisms_resolve_targets_to_symbols(tmp_path):
    """Each drug-target link becomes a gene-drug row named from drug_molecule."""
    mechanisms = [_mechanism(targets=["ENSG01", "ENSG02", "ENSG99"])]
    release = _release(tmp_path, [_report()], mechanisms)

    result = opentargets.parse_mechanisms(
        release["drug_mechanism_of_action_path"],
        release["drug_molecule_path"],
        release["target_path"],
    ).sort_values("gene")

    assert result[["gene", "drug_name", "action_type", "trade_names"]].values.tolist() == [
        ["GENEA", "DRUG1", "INHIBITOR", "Brand1"],
        ["GENEB", "DRUG1", "INHIBITOR", "Brand1"],
    ]


def test_mechanisms_carry_target_class_and_reference_sources(tmp_path):
    """target_class is the top classification level; subsource the mechanism's references."""
    references = [
        {"source": "PubMed", "ids": ["1"], "urls": None},
        {"source": "DailyMed", "ids": ["2"], "urls": None},
        {"source": "PubMed", "ids": ["3"], "urls": None},
    ]
    mechanisms = [_mechanism(targets=["ENSG01", "ENSG02"], references=references)]
    release = _release(tmp_path, [_report()], mechanisms)

    result = opentargets.parse_mechanisms(
        release["drug_mechanism_of_action_path"],
        release["drug_molecule_path"],
        release["target_path"],
    ).sort_values("gene")

    assert result[["gene", "target_class", "subsource"]].values.tolist() == [
        ["GENEA", "Enzyme", "DailyMed|PubMed"],
        ["GENEB", None, "DailyMed|PubMed"],
    ]


def test_datasets_replacing_known_drug_write_the_same_outputs(tmp_path):
    """The four 26.06 datasets write the same two outputs known_drug does."""
    release = tmp_path / "release"
    release.mkdir()

    gene_drug, drug_indication = opentargets.parse_opentargets(
        tmp_path / "out" / "opentargets_gene_drug.parquet",
        tmp_path / "out" / "opentargets_drug_indication.parquet",
        **_release(release, [_report()]),
    )

    assert gene_drug["gene"].tolist() == ["GENEA"]
    assert gene_drug.columns.tolist() == [
        "gene",
        "drug_name",
        "chembl_id",
        "molecule_type",
        "action_type",
        "target_class",
        "trade_names",
        "synonyms",
        "subsource",
    ]
    assert drug_indication[["chembl_id", "efo_id", "phase"]].values.tolist() == [
        ["CHEMBL1", "MONDO:0004979", 2],
    ]
    assert drug_indication["phase"].dtype == "Int64"


@pytest.mark.parametrize(
    "datasets",
    [
        ["clinical_report_path", "drug_molecule_path"],
        [],
        ["known_drug_path", "clinical_report_path"],
    ],
)
def test_known_drug_or_all_four_replacements_are_required(tmp_path, datasets):
    """Exactly one of known_drug or all four of its replacements is needed."""
    with pytest.raises(ValueError, match="known_drug"):
        opentargets.parse_opentargets(
            tmp_path / "gene_drug.parquet",
            tmp_path / "drug_indication.parquet",
            **dict.fromkeys(datasets, tmp_path),
        )
