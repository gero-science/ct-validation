"""Enrichment calculation for clinical trial phase transitions."""

from dataclasses import dataclass

import numpy as np
import pandas as pd

from ct_validation.config.schema import phase_label
from ct_validation.validation.statistics import (
    BootstrapCI,
    cluster_bootstrap_ci,
    fisher_exact_pvalue,
    katz_ci_risk_ratio,
    woolf_ci_odds_ratio,
)


@dataclass
class EnrichmentResult:
    """Result of enrichment calculation for a phase transition."""

    phase_from: int
    phase_to: int
    label: str
    # Counts
    x_yes: int  # Successes with genetic evidence
    n_yes: int  # Total with genetic evidence
    x_no: int  # Successes without genetic evidence
    n_no: int  # Total without genetic evidence
    # Rates
    rate_yes: float  # Success rate with genetic evidence
    rate_no: float  # Success rate without genetic evidence
    # Risk ratio and CI
    rr: float
    rr_ci_lower: float
    rr_ci_upper: float
    # Odds ratio and CI
    or_: float
    or_ci_lower: float
    or_ci_upper: float
    # Fisher's exact test (two-sided)
    p_value: float
    # Cluster bootstrap CIs; None when not requested
    boot: BootstrapCI | None = None

    def to_dict(self) -> dict:
        """Convert to dictionary for DataFrame creation."""
        result = {
            "phase_from": self.phase_from,
            "phase_to": self.phase_to,
            "phase_label": self.label,
            "x_yes": self.x_yes,
            "n_yes": self.n_yes,
            "x_no": self.x_no,
            "n_no": self.n_no,
            "rate_yes": self.rate_yes,
            "rate_no": self.rate_no,
            "rr": self.rr,
            "rr_ci_lower": self.rr_ci_lower,
            "rr_ci_upper": self.rr_ci_upper,
            "or": self.or_,
            "or_ci_lower": self.or_ci_lower,
            "or_ci_upper": self.or_ci_upper,
            "p_value": self.p_value,
        }
        # Omitted when off so the default output schema is unchanged
        if self.boot is not None:
            result["rr_boot_ci_lower"] = self.boot.rr_ci_lower
            result["rr_boot_ci_upper"] = self.boot.rr_ci_upper
            result["or_boot_ci_lower"] = self.boot.or_ci_lower
            result["or_boot_ci_upper"] = self.boot.or_ci_upper
        return result


def calculate_enrichment(
    df: pd.DataFrame,
    phase_from: int,
    phase_to: int,
    phase_col: str = "max_phase",
    genetic_evidence_col: str = "has_genetic_evidence",
    ongoing_col: str = "is_ongoing",
    cluster_col: str | None = None,
    bootstrap_replicates: int = 10_000,
    bootstrap_seed: int = 0,
) -> EnrichmentResult:
    """
    Calculate enrichment for a phase transition.

    The enrichment formula:
        RR = [P(success | genetic_evidence)] / [P(success | no_evidence)]

    Where success = reaching phase_to given starting at phase_from.

    Parameters
    ----------
    df : pd.DataFrame
        Clinical trial data with phase and genetic evidence columns
    phase_from : int
        Starting phase (1-4)
    phase_to : int
        Target phase (must be > phase_from)
    phase_col : str
        Column name for max phase reached
    genetic_evidence_col : str
        Column name for genetic evidence flag
    ongoing_col : str
        Column name for the still-running flag. An absent column censors nothing. Must be
        plain non-null bool: NA censors in opposite directions under object vs
        nullable-boolean dtype. validate() normalizes it, direct callers must do so themselves.
    cluster_col : str | None
        Column to resample as blocks for a secondary bootstrap CI on both the risk and
        the odds ratio, typically "gene". None (default) skips the bootstrap and leaves
        `boot` unset.
    bootstrap_replicates : int
        Bootstrap replicates
    bootstrap_seed : int
        Seed for the bootstrap RNG

    Returns
    -------
    EnrichmentResult
        Enrichment metrics with counts, rates, and confidence intervals
    """
    # Filter to programs that reached at least phase_from
    df_eligible = df[df[phase_col] >= phase_from].copy()

    # A pair still running below phase_to has not failed the transition — its outcome is
    # undetermined, so it belongs in neither numerator nor denominator. Censoring is
    # per-transition: one that already reached phase_to succeeded regardless.
    if ongoing_col in df_eligible.columns:
        undetermined = df_eligible[ongoing_col] & (df_eligible[phase_col] < phase_to)
        df_eligible = df_eligible[~undetermined]

    # Split by genetic evidence
    with_ge = df_eligible[df_eligible[genetic_evidence_col]]
    without_ge = df_eligible[~df_eligible[genetic_evidence_col]]

    # Count successes (reached phase_to)
    n_yes = len(with_ge)
    x_yes = int((with_ge[phase_col] >= phase_to).sum())
    n_no = len(without_ge)
    x_no = int((without_ge[phase_col] >= phase_to).sum())

    # Calculate rates
    rate_yes = x_yes / n_yes if n_yes > 0 else 0.0
    rate_no = x_no / n_no if n_no > 0 else 0.0

    # RR/OR require both comparison groups; continuity correction is not meaningful here
    boot = None
    if n_yes == 0 or n_no == 0:
        rr = rr_lower = rr_upper = np.nan
        or_ = or_lower = or_upper = np.nan
        p_value = np.nan
        if cluster_col is not None:
            boot = BootstrapCI.undefined()
    else:
        rr, rr_lower, rr_upper = katz_ci_risk_ratio(x_yes, n_yes, x_no, n_no)
        or_, or_lower, or_upper = woolf_ci_odds_ratio(x_yes, n_yes, x_no, n_no)
        p_value = fisher_exact_pvalue(x_yes, n_yes, x_no, n_no)

        if cluster_col is not None:
            boot = cluster_bootstrap_ci(
                success=(df_eligible[phase_col] >= phase_to).to_numpy(),
                has_evidence=df_eligible[genetic_evidence_col].to_numpy(dtype=bool),
                cluster=df_eligible[cluster_col].to_numpy(),
                n_replicates=bootstrap_replicates,
                seed=bootstrap_seed,
            )

    # Generate label
    label = phase_label(phase_from, phase_to)

    return EnrichmentResult(
        phase_from=phase_from,
        phase_to=phase_to,
        label=label,
        x_yes=x_yes,
        n_yes=n_yes,
        x_no=x_no,
        n_no=n_no,
        rate_yes=rate_yes,
        rate_no=rate_no,
        rr=rr,
        rr_ci_lower=rr_lower,
        rr_ci_upper=rr_upper,
        or_=or_,
        or_ci_lower=or_lower,
        or_ci_upper=or_upper,
        p_value=p_value,
        boot=boot,
    )


def calculate_all_enrichments(
    df: pd.DataFrame,
    phase_transitions: list[tuple[int, int]],
    phase_col: str = "max_phase",
    genetic_evidence_col: str = "has_genetic_evidence",
    ongoing_col: str = "is_ongoing",
    cluster_col: str | None = None,
    bootstrap_replicates: int = 10_000,
    bootstrap_seed: int = 0,
) -> pd.DataFrame:
    """
    Calculate enrichment for multiple phase transitions.

    Parameters
    ----------
    df : pd.DataFrame
        Clinical trial data
    phase_transitions : list of (from, to) tuples
        Phase transitions to calculate
    phase_col : str
        Column name for max phase
    genetic_evidence_col : str
        Column name for genetic evidence flag
    ongoing_col : str
        Column name for the still-running flag (optional, see calculate_enrichment)
    cluster_col : str | None
        Column to block-resample for a secondary bootstrap CI (see calculate_enrichment)
    bootstrap_replicates : int
        Bootstrap replicates
    bootstrap_seed : int
        Seed for the bootstrap RNG

    Returns
    -------
    pd.DataFrame
        Enrichment results with one row per transition
    """
    results = []
    for phase_from, phase_to in phase_transitions:
        result = calculate_enrichment(
            df,
            phase_from,
            phase_to,
            phase_col=phase_col,
            genetic_evidence_col=genetic_evidence_col,
            ongoing_col=ongoing_col,
            cluster_col=cluster_col,
            bootstrap_replicates=bootstrap_replicates,
            bootstrap_seed=bootstrap_seed,
        )
        results.append(result.to_dict())

    return pd.DataFrame(results)
