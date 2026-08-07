"""Tests for katz_ci_risk_ratio."""

import math

import numpy as np
import pytest
from ct_validation.validation import statistics
from ct_validation.validation.statistics import (
    BootstrapReplicateLossWarning,
    cluster_bootstrap_ci,
    fisher_exact_pvalue,
    katz_ci_risk_ratio,
    woolf_ci_odds_ratio,
)


def test_fisher_exact_equal_rates():
    """p-value is high when both groups have the same success rate."""
    p_value = fisher_exact_pvalue(50, 100, 50, 100)
    assert p_value == pytest.approx(1.0, abs=1e-9)


def test_fisher_exact_enrichment_signal():
    """p-value is low when group 1 has higher success rate."""
    p_value = fisher_exact_pvalue(80, 100, 40, 100)
    assert p_value < 0.05


def test_fisher_exact_empty_group_returns_nan():
    """p-value is undefined when either group is empty."""
    assert math.isnan(fisher_exact_pvalue(0, 0, 50, 100))
    assert math.isnan(fisher_exact_pvalue(10, 10, 0, 0))


def test_rr_equal_rates():
    """RR=1.0 when both groups have same success rate."""
    rr, ci_lower, ci_upper = katz_ci_risk_ratio(50, 100, 50, 100)
    assert rr == pytest.approx(1.0)


def test_rr_higher_in_group1():
    """RR>1.0 when group1 has higher success rate."""
    rr, _, _ = katz_ci_risk_ratio(80, 100, 40, 100)
    assert rr > 1.0
    assert rr == pytest.approx(2.0)


def test_rr_lower_in_group1():
    """RR<1.0 when group1 has lower success rate."""
    rr, _, _ = katz_ci_risk_ratio(20, 100, 40, 100)
    assert rr < 1.0
    assert rr == pytest.approx(0.5)


def test_ci_contains_rr():
    """CI_lower <= RR <= CI_upper."""
    rr, ci_lower, ci_upper = katz_ci_risk_ratio(30, 100, 50, 100)
    assert ci_lower <= rr <= ci_upper


def test_zero_successes_applies_correction():
    """Handles a=0 with continuity correction."""
    rr, ci_lower, ci_upper = katz_ci_risk_ratio(0, 100, 50, 100)
    assert not math.isnan(rr)
    assert rr < 1.0


def test_zero_in_both_groups_applies_correction():
    """Handles a=0 and c=0 with continuity correction."""
    rr, ci_lower, ci_upper = katz_ci_risk_ratio(0, 100, 0, 100)
    assert not math.isnan(rr)
    assert rr == pytest.approx(1.0)


def test_ci_width_decreases_with_sample_size():
    """Larger samples produce narrower CI."""
    _, ci_lower_small, ci_upper_small = katz_ci_risk_ratio(5, 10, 5, 10)
    _, ci_lower_large, ci_upper_large = katz_ci_risk_ratio(500, 1000, 500, 1000)

    width_small = ci_upper_small - ci_lower_small
    width_large = ci_upper_large - ci_lower_large
    assert width_large < width_small


def test_rr_perfect_success_applies_correction():
    """a==n1 (100% success) triggers correction and produces non-degenerate CI."""
    rr, ci_lower, ci_upper = katz_ci_risk_ratio(1, 1, 50, 100)
    assert not math.isnan(rr)
    assert rr > 1.0
    assert ci_lower < rr < ci_upper  # CI must have nonzero width


def test_or_perfect_success_applies_correction():
    """a==n1 (100% success, b=0) triggers correction and produces non-degenerate CI."""
    or_, ci_lower, ci_upper = woolf_ci_odds_ratio(1, 1, 50, 100)
    assert not math.isnan(or_)
    assert ci_lower < or_ < ci_upper


# --- Odds Ratio tests ---


def test_or_equal_rates():
    """OR=1.0 when both groups have same success rate."""
    or_, ci_lower, ci_upper = woolf_ci_odds_ratio(50, 100, 50, 100)
    assert or_ == pytest.approx(1.0)


def test_or_higher_in_group1():
    """OR>1.0 when group1 has higher success rate."""
    or_, _, _ = woolf_ci_odds_ratio(80, 100, 40, 100)
    assert or_ > 1.0
    # OR = (80*60)/(20*40) = 6.0
    assert or_ == pytest.approx(6.0)


def test_or_lower_in_group1():
    """OR<1.0 when group1 has lower success rate."""
    or_, _, _ = woolf_ci_odds_ratio(20, 100, 40, 100)
    assert or_ < 1.0


def test_or_ci_contains_or():
    """CI_lower <= OR <= CI_upper."""
    or_, ci_lower, ci_upper = woolf_ci_odds_ratio(30, 100, 50, 100)
    assert ci_lower <= or_ <= ci_upper


def test_or_zero_successes_applies_correction():
    """Handles a=0 with continuity correction."""
    or_, ci_lower, ci_upper = woolf_ci_odds_ratio(0, 100, 50, 100)
    assert not math.isnan(or_)
    assert or_ < 1.0


def test_or_zero_in_both_groups():
    """Handles a=0 and c=0 with continuity correction."""
    or_, ci_lower, ci_upper = woolf_ci_odds_ratio(0, 100, 0, 100)
    assert not math.isnan(or_)
    assert or_ == pytest.approx(1.0)


def test_or_ci_width_decreases_with_sample_size():
    """Larger samples produce narrower OR CI."""
    _, ci_lower_small, ci_upper_small = woolf_ci_odds_ratio(5, 10, 5, 10)
    _, ci_lower_large, ci_upper_large = woolf_ci_odds_ratio(500, 1000, 500, 1000)

    width_small = ci_upper_small - ci_lower_small
    width_large = ci_upper_large - ci_lower_large
    assert width_large < width_small


def _make_rows(x_yes: int, n_yes: int, x_no: int, n_no: int):
    """Expand 2x2 counts into per-row success/evidence arrays."""
    success = np.array(
        [True] * x_yes + [False] * (n_yes - x_yes) + [True] * x_no + [False] * (n_no - x_no),
    )
    has_evidence = np.array([True] * n_yes + [False] * n_no)
    return success, has_evidence


def _log_width(ci_lower: float, ci_upper: float) -> float:
    """CI width on the log scale, where a Katz interval is symmetric."""
    return math.log(ci_upper) - math.log(ci_lower)


def _rr_log_width(boot) -> float:
    return _log_width(boot.rr_ci_lower, boot.rr_ci_upper)


def test_bootstrap_matches_closed_form_when_every_row_is_its_own_cluster():
    """With no clustering to exploit, the bootstrap reproduces Katz and Woolf."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    cluster = np.arange(len(success))

    _, katz_lower, katz_upper = katz_ci_risk_ratio(200, 500, 100, 500)
    _, woolf_lower, woolf_upper = woolf_ci_odds_ratio(200, 500, 100, 500)
    boot = cluster_bootstrap_ci(success, has_evidence, cluster)

    assert boot.rr_ci_lower == pytest.approx(katz_lower, rel=0.03)
    assert boot.rr_ci_upper == pytest.approx(katz_upper, rel=0.03)
    assert boot.or_ci_lower == pytest.approx(woolf_lower, rel=0.03)
    assert boot.or_ci_upper == pytest.approx(woolf_upper, rel=0.03)


def test_bootstrap_ignores_duplication_within_clusters_while_katz_shrinks():
    """Repeating each cluster's rows adds rows but no independent information."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    cluster = np.arange(len(success))

    reps = 10
    success_dup = np.repeat(success, reps)
    has_evidence_dup = np.repeat(has_evidence, reps)
    cluster_dup = np.repeat(cluster, reps)

    boot = cluster_bootstrap_ci(success, has_evidence, cluster)
    boot_dup = cluster_bootstrap_ci(success_dup, has_evidence_dup, cluster_dup)
    _, katz_lower, katz_upper = katz_ci_risk_ratio(200, 500, 100, 500)
    _, katz_dup_lower, katz_dup_upper = katz_ci_risk_ratio(2000, 5000, 1000, 5000)

    # Katz treats the copies as new samples and narrows by ~sqrt(reps)
    katz_shrink = _log_width(katz_lower, katz_upper) / _log_width(katz_dup_lower, katz_dup_upper)
    assert katz_shrink == pytest.approx(math.sqrt(reps), rel=0.05)

    # Every count in a replicate scales by exactly reps, so the ratios are unchanged
    assert _rr_log_width(boot_dup) == pytest.approx(_rr_log_width(boot), rel=1e-9)


def test_bootstrap_widens_when_outcomes_are_correlated_within_clusters():
    """Clusters that share an outcome carry less information than their row count suggests."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    independent = np.arange(len(success))
    # 10 rows per cluster, each cluster drawn from one arm with one shared outcome
    correlated = independent // 10

    boot_independent = cluster_bootstrap_ci(success, has_evidence, independent)
    boot_correlated = cluster_bootstrap_ci(success, has_evidence, correlated)

    assert _rr_log_width(boot_correlated) > _rr_log_width(boot_independent)
    assert _log_width(boot_correlated.or_ci_lower, boot_correlated.or_ci_upper) > _log_width(
        boot_independent.or_ci_lower,
        boot_independent.or_ci_upper,
    )


def test_bootstrap_is_reproducible_across_seeds():
    """Same seed gives identical bounds; a different seed stays close."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    cluster = np.arange(len(success))

    first = cluster_bootstrap_ci(success, has_evidence, cluster, seed=1)
    repeat = cluster_bootstrap_ci(success, has_evidence, cluster, seed=1)
    other = cluster_bootstrap_ci(success, has_evidence, cluster, seed=2)

    assert first == repeat
    assert other.rr_ci_lower == pytest.approx(first.rr_ci_lower, rel=0.05)
    assert other.rr_ci_upper == pytest.approx(first.rr_ci_upper, rel=0.05)
    assert other.or_ci_lower == pytest.approx(first.or_ci_lower, rel=0.05)
    assert other.or_ci_upper == pytest.approx(first.or_ci_upper, rel=0.05)


def test_bootstrap_accepts_non_boolean_flags():
    """0/1 integer columns resample as flags, not as row indices."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    cluster = np.arange(len(success))

    from_bool = cluster_bootstrap_ci(success, has_evidence, cluster)
    from_int = cluster_bootstrap_ci(success.astype(int), has_evidence.astype(int), cluster)

    assert from_int == from_bool


def test_bootstrap_empty_arm_returns_nan():
    """No comparison group means no ratio to resample."""
    success, has_evidence = _make_rows(200, 500, 0, 0)
    cluster = np.arange(len(success))

    boot = cluster_bootstrap_ci(success, has_evidence, cluster)

    assert math.isnan(boot.rr_ci_lower)
    assert math.isnan(boot.rr_ci_upper)
    assert math.isnan(boot.or_ci_lower)
    assert math.isnan(boot.or_ci_upper)


def test_bootstrap_single_cluster_per_arm_returns_nan():
    """One cluster per arm rebuilds the same table every replicate, so the CI is undefined."""
    success, has_evidence = _make_rows(3, 5, 1, 5)
    cluster = np.where(has_evidence, "GENE_A", "GENE_B")

    boot = cluster_bootstrap_ci(success, has_evidence, cluster)

    assert math.isnan(boot.rr_ci_lower)
    assert math.isnan(boot.rr_ci_upper)
    assert math.isnan(boot.or_ci_lower)
    assert math.isnan(boot.or_ci_upper)


def test_bootstrap_no_rows_returns_nan():
    """Empty input produces no clusters to resample."""
    empty = np.array([], dtype=bool)
    boot = cluster_bootstrap_ci(empty, empty, np.array([]))

    assert math.isnan(boot.rr_ci_lower)
    assert math.isnan(boot.rr_ci_upper)
    assert math.isnan(boot.or_ci_lower)
    assert math.isnan(boot.or_ci_upper)


def test_bootstrap_tolerates_null_cluster_labels():
    """A missing label groups with the other missing labels instead of raising."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    cluster = np.arange(len(success)).astype(object)
    cluster[::10] = None

    boot = cluster_bootstrap_ci(success, has_evidence, cluster)

    assert boot.rr_ci_lower < boot.rr_ci_upper
    assert boot.or_ci_lower < boot.or_ci_upper


def _clustered(per_cluster: int, n_clusters: int, x: int, first_label: int):
    """n_clusters clusters of per_cluster rows, x of them successes."""
    success, cluster = [], []
    for i in range(n_clusters):
        success += [True] * x + [False] * (per_cluster - x)
        cluster += [first_label + i] * per_cluster
    return np.array(success), np.array(cluster)


def test_bootstrap_zero_cell_in_observed_table_returns_nan():
    """A cell that is zero before resampling is zero in every replicate.

    The ratio then varies only with how many clusters the draw put in the other arm, which
    would report an interval narrower than the closed form it exists to widen.
    """
    s_yes, c_yes = _clustered(5, 4, 3, 0)
    s_no, c_no = _clustered(5, 4, 0, 100)  # no control success anywhere
    success = np.concatenate([s_yes, s_no])
    cluster = np.concatenate([c_yes, c_no])
    has_evidence = np.concatenate([np.ones(len(s_yes), bool), np.zeros(len(s_no), bool)])

    boot = cluster_bootstrap_ci(success, has_evidence, cluster)

    assert math.isnan(boot.rr_ci_lower)
    assert math.isnan(boot.or_ci_upper)


def _two_arm(n_evidence: int, x_evidence: int, n_control: int, control_with_success: int):
    """Two arms of two-row clusters; `control_with_success` control clusters carry a success."""
    rows = max(2, x_evidence)
    success, cluster, has_evidence = [], [], []
    for i in range(n_evidence):
        success += [True] * x_evidence + [False] * (rows - x_evidence)
        cluster += [i] * rows
        has_evidence += [True] * rows
    for j in range(n_control):
        x = 1 if j < control_with_success else 0
        success += [True] * x + [False] * (2 - x)
        cluster += [1000 + j] * 2
        has_evidence += [False] * 2
    return np.array(success), np.array(has_evidence), np.array(cluster)


def test_bootstrap_reports_an_infinite_bound_when_a_cell_empties_in_many_replicates():
    """A cell held by one cluster vanishes from a third of draws, sending the ratio to infinity.

    Dropping them would remove only the largest ratios, and report a finite upper bound
    the data does not support.
    """
    success, has_evidence, cluster = _two_arm(12, 1, 12, 1)

    boot = cluster_bootstrap_ci(success, has_evidence, cluster)

    assert math.isinf(boot.rr_ci_upper)
    assert math.isinf(boot.or_ci_upper)
    assert math.isfinite(boot.rr_ci_lower)
    assert math.isfinite(boot.or_ci_lower)


def test_bootstrap_scale_invariance_holds_even_when_replicates_degenerate():
    """Repeating rows within clusters scales every cell alike, so the ratios cannot move.

    A continuity correction would break this: 0.5 added to a resampled cell does not scale
    with it. The property therefore fails for any scheme that corrects replicates, whether
    per replicate or all at once.
    """
    success, has_evidence, cluster = _two_arm(12, 1, 12, 2)
    reps = 7

    plain = cluster_bootstrap_ci(success, has_evidence, cluster)
    scaled = cluster_bootstrap_ci(
        np.repeat(success, reps),
        np.repeat(has_evidence, reps),
        np.repeat(cluster, reps),
    )

    # The fixture degenerates: a cell empties often enough that the upper bound is infinite.
    assert math.isinf(plain.rr_ci_upper) and math.isinf(scaled.rr_ci_upper)
    assert scaled.rr_ci_lower == plain.rr_ci_lower
    assert scaled.or_ci_lower == plain.or_ci_lower


def test_bootstrap_rejects_a_replicate_count_below_one():
    """Zero replicates has no answer to give, and must not fail obscurely downstream."""
    success, has_evidence, cluster = _two_arm(10, 1, 10, 5)

    with pytest.raises(ValueError, match="n_replicates"):
        cluster_bootstrap_ci(success, has_evidence, cluster, n_replicates=0)


def test_bootstrap_warns_when_many_replicates_lose_an_arm():
    """Few clusters means many draws empty an arm, leaving 0/0 and no ratio to contribute."""
    s_yes, c_yes = _clustered(5, 2, 3, 0)
    s_no, c_no = _clustered(5, 8, 1, 100)
    success = np.concatenate([s_yes, s_no])
    cluster = np.concatenate([c_yes, c_no])
    has_evidence = np.concatenate([np.ones(len(s_yes), bool), np.zeros(len(s_no), bool)])

    with pytest.warns(BootstrapReplicateLossWarning, match="undefined"):
        cluster_bootstrap_ci(success, has_evidence, cluster)


def test_bootstrap_result_is_independent_of_draw_chunking(monkeypatch):
    """Chunking bounds peak memory; it must not change the numbers."""
    success, has_evidence = _make_rows(200, 500, 100, 500)
    cluster = np.arange(len(success)) // 4

    whole = cluster_bootstrap_ci(success, has_evidence, cluster, n_replicates=500)
    monkeypatch.setattr(statistics, "_DRAW_CHUNK_BYTES", 8_000)
    chunked = cluster_bootstrap_ci(success, has_evidence, cluster, n_replicates=500)

    assert chunked == whole
