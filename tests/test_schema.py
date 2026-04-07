"""Tests for data schema validation."""

import pandas as pd
import pytest
from ct_validation.data.schema import CLINICAL_TRIALS, GENETIC_EVIDENCE, SIMILARITY_LOOKUP


def test_valid_clinical_trials_passes():
    """DataFrame with required columns passes validation."""
    df = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "efo_id": ["EFO:001"],
            "max_phase": [2],
        }
    )
    # Should not raise
    CLINICAL_TRIALS.validate(df, "test")


def test_valid_genetic_evidence_passes():
    """DataFrame with required columns passes validation."""
    df = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "efo_id": ["EFO:001"],
        }
    )
    GENETIC_EVIDENCE.validate(df, "test")


def test_valid_similarity_lookup_passes():
    """DataFrame with required columns passes validation."""
    df = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001"],
            "efo_id_2": ["EFO:002"],
            "similarity": [0.9],
        }
    )
    SIMILARITY_LOOKUP.validate(df, "test")


def test_missing_column_raises():
    """Missing required column raises ValueError."""
    df = pd.DataFrame(
        {
            "gene": ["GENE1"],
            # missing efo_id and max_phase
        }
    )
    with pytest.raises(ValueError, match="missing required columns"):
        CLINICAL_TRIALS.validate(df, "test")


def test_extra_columns_allowed():
    """Additional columns beyond required are allowed."""
    df = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "efo_id": ["EFO:001"],
            "max_phase": [2],
            "extra_column": ["foo"],
            "another_extra": [123],
        }
    )
    # Should not raise
    CLINICAL_TRIALS.validate(df, "test")


def test_parquet_validation(tmp_path):
    """Validates parquet file metadata."""
    df = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "efo_id": ["EFO:001"],
            "max_phase": [2],
        }
    )
    parquet_path = tmp_path / "test.parquet"
    df.to_parquet(parquet_path)

    CLINICAL_TRIALS.validate_parquet(parquet_path, "test")


def test_parquet_missing_column_raises(tmp_path):
    """Parquet missing required column raises ValueError."""
    df = pd.DataFrame({"gene": ["GENE1"]})  # missing efo_id, max_phase
    parquet_path = tmp_path / "test.parquet"
    df.to_parquet(parquet_path)

    with pytest.raises(ValueError, match="missing required columns"):
        CLINICAL_TRIALS.validate_parquet(parquet_path, "test")
