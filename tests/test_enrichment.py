"""Tests for enrichment calculation."""

import pandas as pd
from ct_validation.validation.enrichment import (
    EnrichmentResult,
    calculate_all_enrichments,
    calculate_enrichment,
)


def test_enrichment_counts_match_input(enrichment_input_df):
    """x_yes, n_yes counts match filtered input."""
    result = calculate_enrichment(enrichment_input_df, phase_from=1, phase_to=2)

    # n_yes = rows with has_genetic_evidence=True and max_phase >= 1
    expected_n_yes = enrichment_input_df[
        (enrichment_input_df["has_genetic_evidence"]) & (enrichment_input_df["max_phase"] >= 1)
    ].shape[0]
    assert result.n_yes == expected_n_yes


def test_enrichment_filters_by_phase_from(enrichment_input_df):
    """Only includes rows with phase >= phase_from."""
    result = calculate_enrichment(enrichment_input_df, phase_from=2, phase_to=3)

    # Should exclude phase 1 rows
    total = result.n_yes + result.n_no
    expected = enrichment_input_df[enrichment_input_df["max_phase"] >= 2].shape[0]
    assert total == expected


def test_enrichment_success_counts_phase_to(enrichment_input_df):
    """x counts rows with phase >= phase_to."""
    result = calculate_enrichment(enrichment_input_df, phase_from=1, phase_to=4)

    # x_yes = GE rows that reached phase 4
    expected_x_yes = enrichment_input_df[
        (enrichment_input_df["has_genetic_evidence"]) & (enrichment_input_df["max_phase"] >= 4)
    ].shape[0]
    assert result.x_yes == expected_x_yes


def test_all_enrichments_returns_dataframe(enrichment_input_df):
    """calculate_all_enrichments returns DataFrame with expected columns."""
    transitions = [(1, 2), (2, 3), (1, 4)]
    result = calculate_all_enrichments(enrichment_input_df, transitions)

    assert isinstance(result, pd.DataFrame)
    assert len(result) == len(transitions)
    assert "phase_label" in result.columns
    assert "rr" in result.columns
    assert "n_yes" in result.columns


def test_no_matching_rows_returns_zero_counts():
    """DataFrame with no rows matching phase_from produces zero counts."""
    df = pd.DataFrame(
        {
            "max_phase": [1, 1, 1],  # No rows >= phase_from=2
            "has_genetic_evidence": [True, False, True],
        }
    )
    result = calculate_enrichment(df, phase_from=2, phase_to=3)

    assert result.n_yes == 0
    assert result.n_no == 0


def test_empty_comparison_group_returns_nan_effect_sizes():
    """RR/OR are undefined when either comparison group is empty."""
    all_with_ge = pd.DataFrame(
        {
            "max_phase": [2, 2, 2],
            "has_genetic_evidence": [True, True, True],
        }
    )
    all_without_ge = pd.DataFrame(
        {
            "max_phase": [2, 2, 2],
            "has_genetic_evidence": [False, False, False],
        }
    )

    with_ge_only = calculate_enrichment(all_with_ge, phase_from=1, phase_to=3)
    without_ge_only = calculate_enrichment(all_without_ge, phase_from=1, phase_to=3)
    no_eligible = calculate_enrichment(all_with_ge, phase_from=3, phase_to=4)

    assert pd.isna(with_ge_only.rr)
    assert pd.isna(with_ge_only.or_)
    assert pd.isna(with_ge_only.p_value)
    assert pd.isna(without_ge_only.rr)
    assert pd.isna(without_ge_only.or_)
    assert pd.isna(without_ge_only.p_value)
    assert pd.isna(no_eligible.rr)
    assert pd.isna(no_eligible.or_)
    assert pd.isna(no_eligible.p_value)


def test_enrichment_includes_p_value(enrichment_input_df):
    """Enrichment results include a Fisher's exact test p-value."""
    result = calculate_enrichment(enrichment_input_df, phase_from=1, phase_to=2)

    assert 0.0 <= result.p_value <= 1.0


def test_all_enrichments_includes_p_value_column(enrichment_input_df):
    """Batch enrichment output includes p_value column."""
    result = calculate_all_enrichments(enrichment_input_df, [(1, 2)])

    assert "p_value" in result.columns


def test_result_to_dict():
    """EnrichmentResult.to_dict() returns expected keys."""
    result = EnrichmentResult(
        phase_from=1,
        phase_to=2,
        label="I→II",
        x_yes=10,
        n_yes=20,
        x_no=30,
        n_no=60,
        rate_yes=0.5,
        rate_no=0.5,
        rr=1.0,
        rr_ci_lower=0.8,
        rr_ci_upper=1.2,
        or_=1.0,
        or_ci_lower=0.7,
        or_ci_upper=1.3,
        p_value=0.42,
    )
    d = result.to_dict()

    assert d["phase_from"] == 1
    assert d["phase_to"] == 2
    assert d["phase_label"] == "I→II"
    assert d["rr"] == 1.0
    assert d["or"] == 1.0
    assert d["or_ci_lower"] == 0.7
    assert d["p_value"] == 0.42


def test_rates_between_zero_and_one(enrichment_input_df):
    """Success rates are in [0, 1]."""
    result = calculate_enrichment(enrichment_input_df, phase_from=1, phase_to=2)

    assert 0 <= result.rate_yes <= 1
    assert 0 <= result.rate_no <= 1
