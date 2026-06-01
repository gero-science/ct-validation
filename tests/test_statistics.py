"""Tests for katz_ci_risk_ratio."""

import math

import pytest
from ct_validation.validation.statistics import (
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
