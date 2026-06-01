"""Tests for CLI."""

import pytest
from ct_validation.cli import cli
from typer.testing import CliRunner


@pytest.fixture
def runner():
    """Typer CLI test runner."""
    return CliRunner()


def test_cli_help(runner):
    """--help works without error."""
    result = runner.invoke(cli, ["--help"])
    assert result.exit_code == 0
    assert "Run clinical trial enrichment validation" in result.output


def test_cli_config_only(runner, tmp_config_yaml, tmp_path):
    """CLI runs with config file only."""
    output_dir = tmp_path / "cli_output"

    result = runner.invoke(
        cli,
        [
            "--config",
            str(tmp_config_yaml),
            "--output-dir",
            str(output_dir),
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / "enrichment_results.parquet").exists()


def test_cli_args_only(runner, tmp_parquet_files, tmp_path):
    """CLI runs with explicit args only."""
    output_dir = tmp_path / "cli_output"

    result = runner.invoke(
        cli,
        [
            "--clinical-trials",
            str(tmp_parquet_files["clinical_trials"]),
            "--targets",
            str(tmp_parquet_files["genetic_evidence"]),
            "--similarity-lookup",
            str(tmp_parquet_files["similarity"]),
            "--output-dir",
            str(output_dir),
            "--similarity-threshold",
            "0.8",
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / "enrichment_results.parquet").exists()


@pytest.mark.parametrize(
    "fmt,ext",
    [
        ("parquet", "parquet"),
        ("tsv", "tsv"),
        ("csv", "csv"),
    ],
)
def test_cli_output_formats(runner, tmp_config_yaml, tmp_path, fmt, ext):
    """CLI supports parquet, tsv, csv output formats."""
    output_dir = tmp_path / "cli_output"

    result = runner.invoke(
        cli,
        [
            "--config",
            str(tmp_config_yaml),
            "--output-dir",
            str(output_dir),
            "--format",
            fmt,
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / f"enrichment_results.{ext}").exists()


def test_cli_save_trials(runner, tmp_config_yaml, tmp_path):
    """--save-trials flag produces annotated_trials file."""
    output_dir = tmp_path / "cli_output"

    result = runner.invoke(
        cli,
        [
            "--config",
            str(tmp_config_yaml),
            "--output-dir",
            str(output_dir),
            "--save-trials",
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / "enrichment_results.parquet").exists()
    assert (output_dir / "annotated_trials.parquet").exists()


def test_cli_save_matched_pairs(runner, tmp_config_yaml, tmp_path):
    """--save-matched-pairs flag produces matched_pairs file."""
    output_dir = tmp_path / "cli_output"

    result = runner.invoke(
        cli,
        [
            "--config",
            str(tmp_config_yaml),
            "--output-dir",
            str(output_dir),
            "--save-matched-pairs",
        ],
    )

    assert result.exit_code == 0
    assert (output_dir / "enrichment_results.parquet").exists()
    assert (output_dir / "matched_pairs.parquet").exists()


def test_cli_prints_summary(runner, tmp_config_yaml, tmp_path):
    """CLI prints enrichment summary to stdout."""
    result = runner.invoke(
        cli,
        [
            "--config",
            str(tmp_config_yaml),
            "--output-dir",
            str(tmp_path),
        ],
    )

    assert result.exit_code == 0
    assert "Enrichment Results" in result.output
    assert "Phase" in result.output
    assert "p-value" in result.output


def test_cli_phase_transitions(runner, tmp_config_yaml, tmp_path):
    """--phase-transitions parses and runs correctly."""
    result = runner.invoke(
        cli,
        [
            "--config",
            str(tmp_config_yaml),
            "--output-dir",
            str(tmp_path),
            "--phase-transitions",
            "1,2;2,3",
        ],
    )

    assert result.exit_code == 0
    assert "I→II" in result.output
    assert "II→III" in result.output
