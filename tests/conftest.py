"""Shared test fixtures."""

import pandas as pd
import pytest


@pytest.fixture
def clinical_trials_df() -> pd.DataFrame:
    """10 clinical trials across phases 1-4."""
    return pd.DataFrame(
        {
            "gene": [
                "GENE1",
                "GENE1",
                "GENE2",
                "GENE2",
                "GENE3",
                "GENE3",
                "GENE4",
                "GENE4",
                "GENE5",
                "GENE5",
            ],
            "efo_id": [
                "EFO:001",
                "EFO:002",
                "EFO:001",
                "EFO:003",
                "EFO:001",
                "EFO:004",
                "EFO:002",
                "EFO:005",
                "EFO:003",
                "EFO:006",
            ],
            "max_phase": [1, 2, 2, 3, 1, 4, 3, 2, 4, 1],
        }
    )


@pytest.fixture
def genetic_evidence_df() -> pd.DataFrame:
    """8 gene-indication associations."""
    return pd.DataFrame(
        {
            "gene": ["GENE1", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7"],
            "efo_id": [
                "EFO:001",
                "EFO:010",
                "EFO:001",
                "EFO:020",
                "EFO:002",
                "EFO:030",
                "EFO:001",
                "EFO:001",
            ],
        }
    )


@pytest.fixture
def similarity_lookup_df() -> pd.DataFrame:
    """EFO similarity pairs (0.5-1.0)."""
    return pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:010", "EFO:020", "EFO:030", "EFO:001", "EFO:002"],
            "efo_id_2": ["EFO:001", "EFO:002", "EFO:004", "EFO:003", "EFO:005", "EFO:002"],
            "similarity": [1.0, 0.85, 0.75, 0.6, 0.55, 1.0],
        }
    )


@pytest.fixture
def enrichment_input_df() -> pd.DataFrame:
    """DataFrame ready for enrichment calculation."""
    return pd.DataFrame(
        {
            "max_phase": [1, 1, 2, 2, 3, 3, 4, 4, 1, 2],
            "has_genetic_evidence": [
                True,
                True,
                True,
                False,
                True,
                False,
                True,
                False,
                False,
                False,
            ],
        }
    )


@pytest.fixture
def tmp_parquet_files(tmp_path, clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """Write fixtures to parquet files and return paths."""
    ct_path = tmp_path / "clinical_trials.parquet"
    ge_path = tmp_path / "genetic_evidence.parquet"
    sim_path = tmp_path / "similarity.parquet"

    clinical_trials_df.to_parquet(ct_path)
    genetic_evidence_df.to_parquet(ge_path)
    similarity_lookup_df.to_parquet(sim_path)

    return {"clinical_trials": ct_path, "genetic_evidence": ge_path, "similarity": sim_path}


@pytest.fixture
def tmp_config_yaml(tmp_path, tmp_parquet_files):
    """Create a config YAML file pointing to fixture parquets."""
    config_path = tmp_path / "config.yaml"
    config_content = f"""
data:
  clinical_trials: "{tmp_parquet_files["clinical_trials"]}"
  genetic_evidence: "{tmp_parquet_files["genetic_evidence"]}"
  efo_similarity_lookup: "{tmp_parquet_files["similarity"]}"

output:
  dir: "{tmp_path / "output"}"
  save_trials: false

thresholds:
  semantic_similarity: 0.8

phase_transitions:
  - [1, 2]
  - [2, 3]
  - [3, 4]
  - [1, 4]
"""
    config_path.write_text(config_content)
    return config_path
