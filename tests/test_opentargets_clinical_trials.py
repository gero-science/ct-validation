"""Tests for clinical_trials/opentargets.py's efo_id derivation.

EFO is a merged ontology: MONDO, HP and Orphanet terms are imported under their own CURIE
prefixes. `parse_drug_indication` keeps each disease id as recorded, converting `_` to `:`
and nothing else.
"""

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
