"""Tests for the p-value cutoff in genetic_evidence/genebass.parse_genebass."""

import pandas as pd
import pytest

PHENOCODE = "1"


@pytest.fixture
def raw_genebass_path(tmp_path):
    """Three associations straddling the default 1e-7 exome-wide cutoff."""
    df = pd.DataFrame(
        {
            "Pvalue": [9e-8, 1e-7, 2e-7],
            "gene_symbol": ["BELOW-THRESHOLD-GENE", "AT-THRESHOLD-GENE", "ABOVE-THRESHOLD-GENE"],
            "phenocode": [PHENOCODE, PHENOCODE, PHENOCODE],
            "coding_description": [None, None, None],
            "description": ["trait", "trait", "trait"],
            "annotation": ["pLoF", "pLoF", "pLoF"],
        }
    )
    path = tmp_path / "genebass.parquet"
    df.to_parquet(path, index=False)
    return path


@pytest.fixture
def manifest_path(tmp_path):
    """Manifest mapping the one phenocode to a single EFO term."""
    df = pd.DataFrame(
        {
            "phenocode": [PHENOCODE],
            "coding_description": [None],
            "trait_efos": ["EFO_0000001"],
            "description": ["trait"],
        }
    )
    path = tmp_path / "manifest.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_association_below_pvalue_threshold_is_kept(
    genebass, raw_genebass_path, manifest_path, tmp_path
):
    """A p-value of 9e-8 clears the exome-wide cutoff and is kept."""
    result = genebass.parse_genebass(raw_genebass_path, manifest_path, tmp_path / "out.parquet")

    assert "BELOW-THRESHOLD-GENE" in set(result["gene"])


def test_association_at_exact_pvalue_threshold_is_dropped(
    genebass, raw_genebass_path, manifest_path, tmp_path
):
    """A p-value of exactly 1e-7 does not clear the cutoff and is dropped.

    genebass_preprocess.py applies the same strict comparison upstream, so a row at
    exactly the threshold is already absent from the shipped parquet.
    """
    result = genebass.parse_genebass(raw_genebass_path, manifest_path, tmp_path / "out.parquet")

    assert "AT-THRESHOLD-GENE" not in set(result["gene"])


def test_association_just_above_pvalue_threshold_is_dropped(
    genebass, raw_genebass_path, manifest_path, tmp_path
):
    """A p-value of 2e-7 does not clear the cutoff and is dropped."""
    result = genebass.parse_genebass(raw_genebass_path, manifest_path, tmp_path / "out.parquet")

    assert "ABOVE-THRESHOLD-GENE" not in set(result["gene"])
