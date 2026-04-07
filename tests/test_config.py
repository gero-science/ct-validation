"""Tests for config loading and schema."""

import pytest
from ct_validation.config import load_config
from ct_validation.config.schema import Config, DataPaths, Output, Thresholds, phase_label


def test_load_config_from_yaml(tmp_config_yaml):
    """Parses valid YAML into Config object."""
    config = load_config(tmp_config_yaml)

    assert isinstance(config, Config)
    assert config.data.clinical_trials is not None
    assert config.thresholds.semantic_similarity == 0.8


def test_config_defaults():
    """Config uses defaults when not specified."""
    config = Config()

    assert config.data.clinical_trials is None
    assert config.output.save_trials is False
    assert len(config.phase_transitions) == 4


def test_thresholds_defaults():
    """Thresholds have expected defaults."""
    thresholds = Thresholds()

    assert thresholds.semantic_similarity == 0.8


def test_data_paths_defaults():
    """DataPaths are None by default."""
    paths = DataPaths()

    assert paths.clinical_trials is None
    assert paths.genetic_evidence is None


def test_output_defaults():
    """Output has expected defaults."""
    output = Output()

    assert output.dir is None
    assert output.save_trials is False


@pytest.mark.parametrize(
    "from_phase,to_phase,expected",
    [
        (1, 2, "I→II"),
        (2, 3, "II→III"),
        (3, 4, "III→Approved"),
        (1, 4, "I→Approved"),
    ],
)
def test_phase_label_standard(from_phase, to_phase, expected):
    """Phase labels match expected format."""
    assert phase_label(from_phase, to_phase) == expected


def test_phase_label_preclinical():
    """Phase 0 produces 'Pre' label."""
    assert phase_label(0, 1) == "Pre→I"


def test_invalid_yaml_raises(tmp_path):
    """Invalid YAML syntax raises error."""
    bad_yaml = tmp_path / "bad.yaml"
    bad_yaml.write_text("data:\n  clinical_trials: [invalid")

    with pytest.raises(Exception):
        load_config(bad_yaml)


def test_partial_yaml_uses_defaults(tmp_path):
    """Missing keys in YAML use defaults."""
    partial_yaml = tmp_path / "partial.yaml"
    partial_yaml.write_text("thresholds:\n  semantic_similarity: 0.9\n")

    config = load_config(partial_yaml)

    assert config.thresholds.semantic_similarity == 0.9
