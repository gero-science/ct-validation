"""Statistical functions for enrichment analysis."""

import numpy as np
from scipy import stats


def katz_ci_risk_ratio(
    a: int,
    n1: int,
    c: int,
    n2: int,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """
    Calculate risk ratio and confidence interval using Katz (log) method.

    Parameters
    ----------
    a : int
        Successes in group 1 (with genetic evidence)
    n1 : int
        Total in group 1
    c : int
        Successes in group 2 (without genetic evidence)
    n2 : int
        Total in group 2
    alpha : float
        Significance level (default 0.05 for 95% CI)

    Returns
    -------
    tuple[float, float, float]
        (risk_ratio, ci_lower, ci_upper)
    """
    b = n1 - a
    d = n2 - c

    # Haldane-Anscombe continuity correction for zero cells
    if a == 0 or b == 0 or c == 0 or d == 0:
        a = a + 0.5
        c = c + 0.5
        n1 = n1 + 1
        n2 = n2 + 1

    p1 = a / n1
    p2 = c / n2

    if p2 == 0:
        return np.nan, np.nan, np.nan

    rr = p1 / p2

    # Katz log method
    log_rr = np.log(rr)
    se_log_rr = np.sqrt((1 / a - 1 / n1) + (1 / c - 1 / n2))
    z = stats.norm.ppf(1 - alpha / 2)

    ci_lower = np.exp(log_rr - z * se_log_rr)
    ci_upper = np.exp(log_rr + z * se_log_rr)

    return rr, ci_lower, ci_upper


def woolf_ci_odds_ratio(
    a: int,
    n1: int,
    c: int,
    n2: int,
    alpha: float = 0.05,
) -> tuple[float, float, float]:
    """
    Calculate odds ratio and confidence interval using Woolf (logit) method.

    Parameters
    ----------
    a : int
        Successes in group 1 (with genetic evidence)
    n1 : int
        Total in group 1
    c : int
        Successes in group 2 (without genetic evidence)
    n2 : int
        Total in group 2
    alpha : float
        Significance level (default 0.05 for 95% CI)

    Returns
    -------
    tuple[float, float, float]
        (odds_ratio, ci_lower, ci_upper)
    """
    b = n1 - a
    d = n2 - c

    # Haldane-Anscombe continuity correction for zero cells
    if a == 0 or b == 0 or c == 0 or d == 0:
        a = a + 0.5
        b = b + 0.5
        c = c + 0.5
        d = d + 0.5

    if a == 0 or b == 0 or c == 0 or d == 0:
        return np.nan, np.nan, np.nan

    odds_ratio = (a * d) / (b * c)

    # Woolf logit method
    log_or = np.log(odds_ratio)
    se_log_or = np.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    z = stats.norm.ppf(1 - alpha / 2)

    ci_lower = np.exp(log_or - z * se_log_or)
    ci_upper = np.exp(log_or + z * se_log_or)

    return odds_ratio, ci_lower, ci_upper
