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


def _censoring_df(rows: list[tuple[int, bool]]) -> pd.DataFrame:
    """Build an enrichment input from (max_phase, is_ongoing) rows, all with evidence."""
    return pd.DataFrame(
        {
            "max_phase": [phase for phase, _ in rows],
            "is_ongoing": [ongoing for _, ongoing in rows],
            "has_genetic_evidence": [True] * len(rows),
        }
    )


def test_ongoing_pair_is_excluded_from_the_transition_it_has_not_finished():
    """A pair still running at phase_from is neither a success nor a failure."""
    df = _censoring_df([(1, True), (1, False)])

    result = calculate_enrichment(df, phase_from=1, phase_to=2)

    assert result.n_yes == 1
    assert result.x_yes == 0


def test_ongoing_pair_above_phase_from_is_censored_for_a_later_target_phase():
    """Censoring keys on phase_to, not phase_from: still running at II is undetermined for IV."""
    df = _censoring_df([(2, True), (2, False)])

    result = calculate_enrichment(df, phase_from=1, phase_to=4)

    assert result.n_yes == 1


def test_ongoing_pair_that_reached_phase_to_still_counts_as_a_success():
    """Reaching the target phase settles the transition however much is still running."""
    df = _censoring_df([(4, True)])

    result = calculate_enrichment(df, phase_from=3, phase_to=4)

    assert result.n_yes == 1
    assert result.x_yes == 1


def test_concluded_pair_below_phase_to_counts_as_a_failure():
    """Only undetermined outcomes are censored; a concluded pair that stalled failed."""
    df = _censoring_df([(2, False)])

    result = calculate_enrichment(df, phase_from=2, phase_to=3)

    assert result.n_yes == 1
    assert result.x_yes == 0


def test_missing_ongoing_column_censors_nothing():
    """Without the column every pair counts as concluded, including the ongoing one."""
    df = _censoring_df([(1, True), (2, False), (4, False)])

    with_col = calculate_enrichment(df, phase_from=1, phase_to=4)
    without_col = calculate_enrichment(df.drop(columns=["is_ongoing"]), phase_from=1, phase_to=4)

    assert with_col.n_yes == 2
    assert without_col.n_yes == 3


def test_censoring_applies_to_every_transition():
    """calculate_all_enrichments passes the flag through to each transition."""
    df = _censoring_df([(1, True), (2, True), (3, False)])

    result = calculate_all_enrichments(df, [(1, 2), (2, 3), (1, 4)])

    n_yes = dict(zip(result["phase_label"], result["n_yes"]))
    assert n_yes["I→II"] == 2  # the phase-1 pair is censored
    assert n_yes["II→III"] == 1  # the phase-2 pair is censored, the phase-1 pair ineligible
    assert n_yes["I→Approved"] == 1  # only the concluded phase-3 pair survives


def _clustered_df() -> pd.DataFrame:
    """40 pairs over 8 genes, with genetic evidence enriched for reaching phase 4.

    Success counts differ between genes; identical clusters would leave the bootstrap no
    variance to report.
    """
    successes_per_gene = {True: [4, 3, 3, 2], False: [2, 1, 1, 0]}
    rows = []
    for has_ge, counts in successes_per_gene.items():
        for gene_idx, n_success in enumerate(counts):
            for pair_idx in range(5):
                rows.append(
                    {
                        "gene": f"GENE_{has_ge}_{gene_idx}",
                        "max_phase": 4 if pair_idx < n_success else 1,
                        "has_genetic_evidence": has_ge,
                    },
                )
    return pd.DataFrame(rows)


def test_bootstrap_columns_absent_by_default():
    """Without cluster_col the output schema is unchanged."""
    result = calculate_all_enrichments(_clustered_df(), [(1, 4)])

    assert "rr_boot_ci_lower" not in result.columns
    assert "rr_boot_ci_upper" not in result.columns


def test_cluster_col_adds_bootstrap_ci_around_both_ratios():
    """cluster_col produces bootstrap intervals bracketing both point estimates."""
    result = calculate_all_enrichments(_clustered_df(), [(1, 4)], cluster_col="gene")

    row = result.iloc[0]
    assert row["rr_boot_ci_lower"] < row["rr"] < row["rr_boot_ci_upper"]
    assert row["or_boot_ci_lower"] < row["or"] < row["or_boot_ci_upper"]


def test_bootstrap_leaves_point_estimates_untouched():
    """Resampling reports uncertainty only; RR, OR and p-value are unaffected."""
    df = _clustered_df()
    plain = calculate_all_enrichments(df, [(1, 4)]).iloc[0]
    clustered = calculate_all_enrichments(df, [(1, 4)], cluster_col="gene").iloc[0]

    for col in ("rr", "rr_ci_lower", "rr_ci_upper", "or", "or_ci_lower", "or_ci_upper", "p_value"):
        assert clustered[col] == plain[col]
