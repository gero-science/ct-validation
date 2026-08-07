"""Statistical functions for enrichment analysis."""

import warnings
from typing import NamedTuple

import numpy as np
import pandas as pd
from scipy import stats


class BootstrapReplicateLossWarning(UserWarning):
    """Raised when many replicates left a ratio undefined (0/0) and were dropped.

    The interval then conditions on the rest. Own class so a caller working with few
    clusters can silence exactly this.
    """


class BootstrapCI(NamedTuple):
    """Bootstrap bounds for both effect measures, from one set of replicates."""

    rr_ci_lower: float
    rr_ci_upper: float
    or_ci_lower: float
    or_ci_upper: float

    @classmethod
    def undefined(cls) -> "BootstrapCI":
        """Bounds for a table resampling cannot speak to."""
        return cls(np.nan, np.nan, np.nan, np.nan)


# Cap on the (n_replicates, n_clusters) int64 draw matrix, which is otherwise the peak
# allocation: the default 10,000 replicates over a 20k-gene universe is 1.6 GB in one block.
_DRAW_CHUNK_BYTES = 64_000_000

# Fraction of replicates with no ratio at all (0/0). Below the first the percentiles barely
# move; above it the caller should know the bounds condition on the rest; past the second
# there is too little left to read an interval off.
_REPLICATE_LOSS_WARN_AT = 0.05
_REPLICATE_LOSS_REFUSE_AT = 0.25


def fisher_exact_pvalue(
    a: int,
    n1: int,
    c: int,
    n2: int,
) -> float:
    """
    Two-sided Fisher's exact test p-value for a 2x2 contingency table.

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

    Returns
    -------
    float
        Two-sided p-value, or NaN when either group is empty
    """
    if n1 == 0 or n2 == 0:
        return np.nan

    b = n1 - a
    d = n2 - c
    _, p_value = stats.fisher_exact([[a, b], [c, d]], alternative="two-sided")
    return float(p_value)


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


def cluster_bootstrap_ci(
    success: np.ndarray,
    has_evidence: np.ndarray,
    cluster: np.ndarray,
    alpha: float = 0.05,
    n_replicates: int = 10_000,
    seed: int = 0,
) -> BootstrapCI:
    """
    Percentile bootstrap CIs for the risk and odds ratios, resampling clusters not rows.

    Katz and Woolf assume independent rows and run too narrow when outcomes correlate
    within a cluster — one gene across many indications. Only the intervals are returned;
    resampling leaves the point estimates unchanged.

    Parameters
    ----------
    success : np.ndarray
        Boolean per row: reached the target phase
    has_evidence : np.ndarray
        Boolean per row: has genetic evidence
    cluster : np.ndarray
        Cluster label per row (e.g. gene)
    alpha : float
        Significance level (default 0.05 for 95% CI)
    n_replicates : int
        Bootstrap replicates
    seed : int
        Seed for the resampling RNG

    Returns
    -------
    BootstrapCI
        Bounds for both measures, all NaN when the table cannot support resampling:
        fewer than two clusters in an arm, a zero cell in the observed table, or too
        many replicates leaving the ratio undefined. A bound may be infinite when the
        resample genuinely cannot bound that side.
    """
    if n_replicates < 1:
        raise ValueError(f"n_replicates must be at least 1, got {n_replicates}")

    success = np.asarray(success, dtype=bool).astype(float)
    has_evidence = np.asarray(has_evidence, dtype=bool)

    # sort=True keys the cluster slots to the labels themselves rather than to the order the
    # rows arrived in, so the draws below are a function of the data alone: the same table
    # sorted differently returns the same interval from the same seed.
    # use_na_sentinel=False so nulls form their own cluster; the default -1 sentinel would
    # make the np.bincount calls below raise. Note that this is a modelling choice and not
    # only a mechanical one: every null-labelled row becomes a single large cluster that
    # dominates the resample, so pass a column without nulls where that matters.
    codes, uniques = pd.factorize(cluster, sort=True, use_na_sentinel=False)
    n_clusters = len(uniques)

    # A one-cluster arm contributes no variance of its own — every replicate rebuilds it at
    # the same rate, since drawing it k times scales a and n1 alike and k cancels. The
    # interval would then be a reading of the other arm only.
    min_clusters_per_arm = 2
    if min(len(np.unique(codes[has_evidence])), len(np.unique(codes[~has_evidence]))) < (
        min_clusters_per_arm
    ):
        return BootstrapCI.undefined()

    a = np.bincount(codes[has_evidence], weights=success[has_evidence], minlength=n_clusters)
    n1 = np.bincount(codes[has_evidence], minlength=n_clusters)
    c = np.bincount(codes[~has_evidence], weights=success[~has_evidence], minlength=n_clusters)
    n2 = np.bincount(codes[~has_evidence], minlength=n_clusters)

    # A zero cell in the observed table survives every resample, so the ratio degenerates
    # into a function of how many clusters the draw happened to put in the other arm. The
    # resulting interval measures the draw rather than the data, and comes out narrower
    # than the closed form it exists to widen — so refuse instead of reporting it.
    if 0 in (a.sum(), n1.sum() - a.sum(), c.sum(), n2.sum() - c.sum()):
        return BootstrapCI.undefined()

    # Drawing n_clusters clusters with replacement is a multinomial over cluster slots, so
    # one matmul per cell replaces a per-replicate resampling loop. Chunked over replicates
    # because the draw matrix, not the results, is what dominates memory.
    rng = np.random.default_rng(seed)
    chunk = max(1, min(n_replicates, _DRAW_CHUNK_BYTES // (n_clusters * 8)))
    parts = []
    for start in range(0, n_replicates, chunk):
        draw_counts = rng.multinomial(
            n_clusters,
            np.full(n_clusters, 1.0 / n_clusters),
            size=min(chunk, n_replicates - start),
        )
        parts.append((draw_counts @ a, draw_counts @ n1, draw_counts @ c, draw_counts @ n2))
    a_b, n1_b, c_b, n2_b = (np.concatenate(x) for x in zip(*parts))

    # No continuity correction here, unlike katz_ci_risk_ratio and woolf_ci_odds_ratio:
    # Haldane-Anscombe exists to keep log(ratio) and its standard error defined, and a
    # percentile bootstrap takes neither. It needs only a ratio per replicate, and a replicate
    # that empties a cell still has one — zero or infinity, order statistics like any other.
    # Discarding those would clip both tails the interval exists to measure: an emptied
    # control cell sends the ratio to infinity, an emptied evidence cell to zero.
    b_b, d_b = n1_b - a_b, n2_b - c_b
    with np.errstate(divide="ignore", invalid="ignore"):
        rr_b = (a_b / n1_b) / (c_b / n2_b)
        or_b = (a_b * d_b) / (b_b * c_b)

    # 0/0 is the one draw with no ratio at all: an arm took no cluster, neither arm kept a
    # success, or — for the odds ratio alone — neither arm kept a failure. Hence a mask per
    # measure: that last case leaves the risk ratio perfectly well defined.
    rr_defined = ~np.isnan(rr_b)
    or_defined = ~np.isnan(or_b)

    lost = 1 - min(rr_defined.mean(), or_defined.mean())
    if lost > _REPLICATE_LOSS_REFUSE_AT:
        return BootstrapCI.undefined()
    if lost > _REPLICATE_LOSS_WARN_AT:
        warnings.warn(
            f"{lost:.1%} of {n_replicates:,} bootstrap replicates left a ratio undefined "
            f"(0/0) and were dropped. The bounds below condition on the rest. More "
            f"clusters, or a table whose cells are spread over more of them, would fix it.",
            BootstrapReplicateLossWarning,
            stacklevel=3,
        )

    # method="lower"/"higher" rather than the default linear interpolation, which returns
    # nan when it has to interpolate onto an infinite neighbour (0 * inf).
    lo_pct, hi_pct = 100 * alpha / 2, 100 * (1 - alpha / 2)
    bounds = []
    for values, defined in ((rr_b, rr_defined), (or_b, or_defined)):
        kept = values[defined]
        bounds += [
            float(np.percentile(kept, lo_pct, method="lower")),
            float(np.percentile(kept, hi_pct, method="higher")),
        ]
    return BootstrapCI(*bounds)


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
