"""Tests for the p-value cutoff in genetic_evidence/gwas_catalog.parse_gwas_catalog."""

import pandas as pd
import pytest

EFO_URI = "http://www.ebi.ac.uk/efo/EFO_0000001"


@pytest.fixture
def raw_gwas_path(tmp_path):
    """Three associations straddling the default 5e-8 genome-wide significance cutoff."""
    df = pd.DataFrame(
        {
            "P-VALUE": [4e-8, 5e-8, 6e-8],
            "MAPPED_GENE": ["BELOW-THRESHOLD-GENE", "AT-THRESHOLD-GENE", "ABOVE-THRESHOLD-GENE"],
            "MAPPED_TRAIT": ["trait", "trait", "trait"],
            "MAPPED_TRAIT_URI": [EFO_URI, EFO_URI, EFO_URI],
        }
    )
    path = tmp_path / "gwas_catalog.tsv"
    df.to_csv(path, sep="\t", index=False)
    return path


def test_association_below_pvalue_threshold_is_kept(gwas_catalog, raw_gwas_path, tmp_path):
    """A p-value of 4e-8 clears genome-wide significance and is kept."""
    result = gwas_catalog.parse_gwas_catalog(raw_gwas_path, tmp_path / "out.parquet")

    assert "BELOW-THRESHOLD-GENE" in set(result["gene"])


def test_association_at_exact_pvalue_threshold_is_dropped(gwas_catalog, raw_gwas_path, tmp_path):
    """A p-value of exactly 5e-8 does not clear genome-wide significance and is dropped."""
    result = gwas_catalog.parse_gwas_catalog(raw_gwas_path, tmp_path / "out.parquet")

    assert "AT-THRESHOLD-GENE" not in set(result["gene"])


def test_association_just_above_pvalue_threshold_is_dropped(gwas_catalog, raw_gwas_path, tmp_path):
    """A p-value of 6e-8 does not clear genome-wide significance and is dropped."""
    result = gwas_catalog.parse_gwas_catalog(raw_gwas_path, tmp_path / "out.parquet")

    assert "ABOVE-THRESHOLD-GENE" not in set(result["gene"])
