"""Tests for EFO loading and obsolete-id remapping."""

import pandas as pd
import pytest

from ct_validation.efo import efo_for_remap, load_efo, remap_efo_ids, to_curie

OBO = """\
format-version: 1.2
data-version: http://www.ebi.ac.uk/efo/releases/v3.90.0/efo.owl

[Term]
id: efo:EFO_0000001
name: live efo

[Term]
id: MONDO:0000001
name: live mondo

[Term]
id: Orphanet:1
name: live orphanet

[Term]
id: efo:EFO_0000010
name: obsolete to mondo by purl
is_obsolete: true
replaced_by: http://purl.obolibrary.org/obo/MONDO_0000001

[Term]
id: efo:EFO_0000011
name: obsolete to efo by iri
is_obsolete: true
replaced_by: http://www.ebi.ac.uk/efo/EFO_0000001

[Term]
id: efo:EFO_0000012
name: obsolete to orphanet by iri
is_obsolete: true
replaced_by: http://www.orpha.net/ORDO/Orphanet_1

[Term]
id: efo:EFO_0000013
name: start of a chain
is_obsolete: true
replaced_by: efo:EFO_0000010

[Term]
id: efo:EFO_0000020
name: consider only
is_obsolete: true
consider: MONDO:0000001

[Term]
id: efo:EFO_0000021
name: no replacement
is_obsolete: true

[Term]
id: efo:EFO_0000022
name: split
is_obsolete: true
replaced_by: MONDO:0000001
replaced_by: Orphanet:1

[Term]
id: efo:EFO_0000023
name: cycle a
is_obsolete: true
replaced_by: efo:EFO_0000024

[Term]
id: efo:EFO_0000024
name: cycle b
is_obsolete: true
replaced_by: efo:EFO_0000023

[Term]
id: efo:EFO_0000025
name: chain to a dead end
is_obsolete: true
replaced_by: efo:EFO_0000021

[Term]
id: efo:EFO_0000026
name: split with one live successor
is_obsolete: true
replaced_by: MONDO:0000001
replaced_by: MONDO:9999999

[Term]
id: efo:EFO_0000027
name: one successor spelled twice
is_obsolete: true
replaced_by: http://purl.obolibrary.org/obo/MONDO_0000001
replaced_by: MONDO:0000001

[Term]
id: efo:EFO_0000028
name: chain through a split with one live successor
is_obsolete: true
replaced_by: efo:EFO_0000026

[Typedef]
id: efo:EFO_9999999
name: a relation, not a term
"""


@pytest.fixture
def efo(tmp_path):
    path = tmp_path / "efo.obo"
    path.write_text(OBO)
    return load_efo(path)


@pytest.mark.parametrize(
    ("term", "expected"),
    [
        ("efo:EFO_0000094", "EFO:0000094"),
        ("http://www.ebi.ac.uk/efo/EFO_0000094", "EFO:0000094"),
        ("http://purl.obolibrary.org/obo/MONDO_0000001", "MONDO:0000001"),
        ("http://www.orpha.net/ORDO/Orphanet_1", "Orphanet:1"),
        ("MONDO:0000001", "MONDO:0000001"),
    ],
)
def test_to_curie(term, expected):
    assert to_curie(term) == expected


def test_live_terms_and_version(efo):
    assert efo.version == "3.90.0"
    assert efo.live == {"EFO:0000001", "MONDO:0000001", "Orphanet:1"}
    assert efo.labels["EFO:0000001"] == "live efo"


def test_replacements_follow_chains(efo):
    assert efo.replacements == {
        "EFO:0000010": "MONDO:0000001",
        "EFO:0000011": "EFO:0000001",
        "EFO:0000012": "Orphanet:1",
        "EFO:0000013": "MONDO:0000001",
        "EFO:0000026": "MONDO:0000001",
        "EFO:0000027": "MONDO:0000001",
        "EFO:0000028": "MONDO:0000001",
    }


def test_unresolved(efo):
    """`consider` only, no replacement, a split into live terms, a cycle and a dead end."""
    assert efo.unresolved == {
        "EFO:0000020",
        "EFO:0000021",
        "EFO:0000022",
        "EFO:0000023",
        "EFO:0000024",
        "EFO:0000025",
    }


def test_remap_efo_ids(efo):
    df = pd.DataFrame(
        {
            "gene": ["A", "B", "C", "D", "E"],
            "efo_id": ["EFO:0000001", "EFO:0000013", "EFO:0000020", "HP:0000001", None],
            "efo_term": ["live efo", "old name", "old name", "outside efo", None],
        },
    )
    out = remap_efo_ids(df, efo, label_column="efo_term")
    assert out["gene"].tolist() == ["A", "B", "D", "E"]
    assert out["efo_id"].tolist() == ["EFO:0000001", "MONDO:0000001", "HP:0000001", None]
    assert out["efo_term"].tolist() == ["live efo", "live mondo", "outside efo", None]
    # The input is not modified
    assert df["efo_id"].iloc[1] == "EFO:0000013"


def test_efo_for_remap(tmp_path):
    path = tmp_path / "efo.obo"
    path.write_text(OBO)
    assert efo_for_remap({"efo_ontology": str(path)}) is None
    assert efo_for_remap({"efo_ontology": str(path), "efo_remap_obsolete": False}) is None
    assert (
        efo_for_remap({"efo_ontology": str(path), "efo_remap_obsolete": True}).version == "3.90.0"
    )
