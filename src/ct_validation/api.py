"""Main API for clinical trial validation."""

import warnings
from dataclasses import dataclass
from pathlib import Path

import duckdb
import pandas as pd

from ct_validation.config import load_config
from ct_validation.config.schema import DEFAULT_PHASE_TRANSITIONS, Config, Thresholds
from ct_validation.data.schema import CLINICAL_TRIALS, GENETIC_EVIDENCE, SIMILARITY_LOOKUP
from ct_validation.validation.enrichment import calculate_all_enrichments
from ct_validation.validation.matching import (
    check_similarity_lookup,
    create_matched_pairs_df,
    create_matched_pairs_set,
)

_DataInput = pd.DataFrame | str | Path

# Defaults from config schema (single source of truth)
_DEFAULTS = Thresholds()


class MissingOngoingWarning(UserWarning):
    """Raised when clinical_trials carries no is_ongoing column.

    Running uncensored is a legitimate choice — it is the convention of Nelson et al. 2015
    and Tsepilov et al. 2026, and the only reading available when the source carries no
    trial status. Its own class so those callers can silence exactly this, rather than
    every UserWarning their program emits.
    """


@dataclass
class _ResolvedParams:
    """Resolved parameters from config and/or explicit args."""

    clinical_trials: _DataInput
    targets_list: list[_DataInput]
    similarity_lookup: _DataInput | None
    baseline_evidence: _DataInput | None
    gene_universe: set[str] | str | Path | None
    similarity_threshold: float
    phase_transitions: list[tuple[int, int]]


def _resolve_params(  # noqa: C901, PLR0912
    config: Config | str | Path | None,
    clinical_trials: _DataInput | None,
    targets: _DataInput | list[_DataInput] | None,
    similarity_lookup: _DataInput | None,
    baseline_evidence: _DataInput | None,
    gene_universe: set[str] | str | Path | None,
    similarity_threshold: float | None,
    phase_transitions: list[tuple[int, int]] | None,
) -> _ResolvedParams:
    """Resolve parameters from config and explicit args (args override config)."""
    # Load config if provided
    if config is not None:
        if isinstance(config, (str, Path)):
            config = load_config(config)

        if clinical_trials is None:
            clinical_trials = config.data.clinical_trials
        if targets is None:
            targets = config.data.genetic_evidence
        if similarity_lookup is None:
            similarity_lookup = config.data.efo_similarity_lookup
        if baseline_evidence is None:
            baseline_evidence = config.data.baseline_evidence

        if gene_universe is None:
            gene_universe = config.data.gene_universe
        if similarity_threshold is None:
            similarity_threshold = config.thresholds.semantic_similarity
        if phase_transitions is None:
            phase_transitions = [tuple(t) for t in config.phase_transitions]

    # Apply defaults from Thresholds (single source of truth)
    if similarity_threshold is None:
        similarity_threshold = _DEFAULTS.semantic_similarity
    if phase_transitions is None:
        phase_transitions = list(DEFAULT_PHASE_TRANSITIONS)

    # Validate required
    if clinical_trials is None:
        raise ValueError("clinical_trials is required")
    if targets is None:
        raise ValueError("targets is required")

    # Normalize targets to list
    targets_list = targets if isinstance(targets, list) else [targets]

    return _ResolvedParams(
        clinical_trials=clinical_trials,
        targets_list=targets_list,
        similarity_lookup=similarity_lookup,
        baseline_evidence=baseline_evidence,
        gene_universe=gene_universe,
        similarity_threshold=similarity_threshold,
        phase_transitions=phase_transitions,
    )


def _validate_input(data: pd.DataFrame | str | Path, schema, name: str) -> None:
    """Validate schema from DataFrame or parquet metadata."""
    if isinstance(data, pd.DataFrame):
        schema.validate(data, name)
    else:
        schema.validate_parquet(data, name)


def _load_gene_universe(source: set[str] | str | Path | None) -> set[str] | None:
    """Load gene universe from set or file."""
    if source is None:
        return None
    if isinstance(source, set):
        return source
    with Path(source).open() as f:
        return {line.strip() for line in f if line.strip()}


def _load_clinical_trials(clinical_trials: pd.DataFrame | str | Path) -> pd.DataFrame:
    """
    Load clinical trials once and enforce the input contract.

    Kept out of the per-target path so batch mode does not re-read and re-warn for
    every target set.
    """
    if isinstance(clinical_trials, pd.DataFrame):
        ct = clinical_trials.copy()
    else:
        ct = pd.read_parquet(clinical_trials)

    ct = ct[ct["gene"].notna() & ct["efo_id"].notna() & ct["max_phase"].notna()]

    # Rows are counted individually, so a repeated pair inflates n. Refuse rather than
    # collapse: only the caller knows what max_phase and is_ongoing mean for a merged row.
    duplicated = ct.duplicated(subset=["gene", "efo_id"])
    if duplicated.any():
        examples = ", ".join(
            f"{gene}/{efo_id}"
            for gene, efo_id in ct.loc[duplicated, ["gene", "efo_id"]].head(3).to_numpy()
        )
        raise ValueError(
            f"clinical_trials has {int(duplicated.sum()):,} duplicate (gene, efo_id) rows "
            f"(e.g. {examples}); aggregate to one row per pair before validating",
        )

    if "is_ongoing" not in ct.columns:
        warnings.warn(
            "clinical_trials has no 'is_ongoing' column, so nothing is censored and a pair "
            "still in progress counts as a failure. Silence this with "
            "warnings.filterwarnings('ignore', category=MissingOngoingWarning).",
            MissingOngoingWarning,
            stacklevel=3,
        )
        return ct

    # A missing flag is not evidence of being ongoing, so it falls back to the missing-column
    # default; otherwise object-dtype None and nullable pd.NA censor in opposite directions.
    offenders = set(ct["is_ongoing"].dropna().unique()) - {True, False}
    if offenders:
        shown = ", ".join(sorted(map(repr, offenders))[:3])
        raise ValueError(
            f"clinical_trials['is_ongoing'] must be boolean (or 0/1), got {shown}",
        )

    # astype before fillna: fillna on an object column downcasts with a FutureWarning
    # under pandas 2.
    return ct.assign(
        is_ongoing=ct["is_ongoing"].astype("boolean").fillna(value=False).astype(bool),
    )


def _prepare_clinical_trials(
    ct: pd.DataFrame,
    target_pairs: set[tuple[str, str]],
    baseline_pairs: set[tuple[str, str]],
    gene_universe: set[str] | None = None,
) -> pd.DataFrame:
    """
    Annotate loaded clinical trials with genetic evidence for enrichment.

    Returns DataFrame with `has_genetic_evidence` column:
    - True = YES group (has matching targets)
    - False = NO group (no targets AND no baseline)

    Rows with baseline-only support are filtered out.
    """
    if gene_universe is not None:
        ct = ct[ct["gene"].isin(gene_universe)]

    pairs = list(zip(ct["gene"], ct["efo_id"]))
    has_target = [p in target_pairs for p in pairs]
    has_baseline = [p in baseline_pairs for p in pairs]

    ct = ct.assign(has_genetic_evidence=has_target, _has_baseline=has_baseline)
    ct = ct[ct["has_genetic_evidence"] | ~ct["_has_baseline"]]

    return ct.drop(columns=["_has_baseline"])


def _pack_validation_result(
    enrichment: pd.DataFrame,
    trials: pd.DataFrame | None,
    matched_pairs: pd.DataFrame | None,
    *,
    return_trials: bool,
    return_matched_pairs: bool,
) -> pd.DataFrame | tuple:
    """Build validate() return value from optional extras."""
    if not return_trials and not return_matched_pairs:
        return enrichment
    parts: list[pd.DataFrame] = [enrichment]
    if return_trials:
        parts.append(trials)
    if return_matched_pairs:
        parts.append(matched_pairs)
    return tuple(parts)


def validate(
    config: Config | str | Path | None = None,
    *,
    clinical_trials: _DataInput | None = None,
    targets: _DataInput | list[_DataInput] | None = None,
    similarity_lookup: _DataInput | None = None,
    baseline_evidence: _DataInput | None = None,
    gene_universe: set[str] | str | Path | None = None,
    similarity_threshold: float | None = None,
    phase_transitions: list[tuple[int, int]] | None = None,
    cluster_col: str | None = None,
    bootstrap_replicates: int = 10_000,
    bootstrap_seed: int = 0,
    return_trials: bool = False,
    return_matched_pairs: bool = False,
    check_similarity: bool = False,
) -> pd.DataFrame | tuple | list:
    """
    Run validation pipeline.

    Can be called with:
    - Config only: validate(config="config.yaml")
    - Args only: validate(clinical_trials=..., targets=..., similarity_lookup=...)
    - Mixed (args override config): validate(config="config.yaml", targets=custom_df)
    - Batch: validate(targets=[gwas_df, clinvar_df, ...], ...) — returns list of results

    Two modes:
    1. **Baseline mode** (baseline_evidence=None):
       - YES: clinical trials with matching targets
       - NO: clinical trials without matching targets

    2. **Prioritized mode** (baseline_evidence provided):
       - YES: clinical trials with matching targets
       - NO: clinical trials without targets AND without baseline_evidence

    Input Schemas:
        clinical_trials: gene, efo_id, max_phase [, is_ongoing]
            Unique gene-indication pairs with max phase reached. is_ongoing censors a pair
            still running below phase_to from that transition; without the column nothing
            is censored and a MissingOngoingWarning is issued.
        targets: gene, efo_id
            Gene-indication pairs with supporting genetic evidence.
            Pass a list of DataFrames/paths for batch validation.
        similarity_lookup: efo_id_1, efo_id_2, similarity (optional)
            Symmetric EFO semantic similarity pairs (diagonal has similarity=1.0).
            Matching probes a single orientation, so the table MUST be symmetric; pass
            check_similarity=True to verify (or call check_similarity_lookup() once
            yourself before a validation loop). If None, performs exact matching on
            (gene, efo_id) without expansion.
        baseline_evidence: gene, efo_id (optional)
            Gene-indication pairs for baseline comparison in prioritized mode.
        check_similarity: bool (default False)
            When True and a similarity_lookup is given, assert it is symmetric (raising
            if not) before matching. Off by default because the check scans the whole
            lookup — for repeated calls, validate once with check_similarity_lookup()
            and leave this False.

    Output Schema:
        phase_from, phase_to: Phase transition (e.g., 1→2, 2→3).
        n_yes, n_no: Pairs at phase_from+ less those censored as still running below
            phase_to (with/without genetic evidence); n is transition-specific.
        x_yes, x_no: Pairs reaching phase_to+ (with/without genetic evidence).
        rate_yes, rate_no: Progression rates.
        rr, rr_ci_lower, rr_ci_upper: Risk ratio with 95% CI (Katz log method).
        or, or_ci_lower, or_ci_upper: Odds ratio with 95% CI (Woolf logit method).
        p_value: Two-sided Fisher's exact test p-value.
        rr_boot_ci_lower/upper, or_boot_ci_lower/upper: only when cluster_col is set —
            percentile bootstrap intervals resampling clusters rather than rows, leaving
            the point estimates unchanged.

    Bootstrap arguments (all ignored unless cluster_col is set):
        cluster_col: Column of clinical_trials to resample as blocks (typically "gene")
            for a secondary bootstrap CI on both ratios, since Katz and Woolf assume
            independent rows and run narrow when one gene spans many indications.
            None (default) skips the bootstrap.
        bootstrap_replicates: Resampling draws behind those intervals (default 10,000).
            The percentile bounds carry Monte Carlo error of their own, falling as
            1/sqrt(n); lower this only for exploratory runs. Ignored without cluster_col.
        bootstrap_seed: Seed for the resampling RNG (default 0). Vary it to check the
            bounds are stable. Ignored without cluster_col.

    Returns:
        Single targets: enrichment DataFrame, or tuple with optional extras in order:
        (enrichment), (enrichment, trials), (enrichment, matched_pairs),
        or (enrichment, trials, matched_pairs).
        List of targets: list of the above, one per target set.
    """
    batch = isinstance(targets, list)

    # Resolve all parameters (normalizes targets to list internally)
    p = _resolve_params(
        config=config,
        clinical_trials=clinical_trials,
        targets=targets,
        similarity_lookup=similarity_lookup,
        baseline_evidence=baseline_evidence,
        gene_universe=gene_universe,
        similarity_threshold=similarity_threshold,
        phase_transitions=phase_transitions,
    )

    # Validate schemas
    _validate_input(p.clinical_trials, CLINICAL_TRIALS, "clinical_trials")
    for t in p.targets_list:
        _validate_input(t, GENETIC_EVIDENCE, "targets")
    if p.similarity_lookup is not None:
        _validate_input(p.similarity_lookup, SIMILARITY_LOOKUP, "similarity_lookup")
        if check_similarity:
            check_similarity_lookup(p.similarity_lookup)
    if p.baseline_evidence is not None:
        _validate_input(p.baseline_evidence, GENETIC_EVIDENCE, "baseline_evidence")

    gu = _load_gene_universe(p.gene_universe)
    ct_all = _load_clinical_trials(p.clinical_trials)

    # Checked here rather than at the bootstrap call: that runs once per phase transition,
    # after matching and pair construction, so a typo would cost the whole pipeline first.
    if cluster_col is not None and cluster_col not in ct_all.columns:
        raise KeyError(
            f"cluster_col={cluster_col!r} is not a column of clinical_trials "
            f"(has {list(ct_all.columns)})",
        )

    # Shared connection for batch efficiency
    con = duckdb.connect()

    # Match baseline once (shared across all target sets)
    baseline_pairs: set[tuple[str, str]] = set()
    if p.baseline_evidence is not None:
        baseline_pairs = create_matched_pairs_set(
            genetic_evidence=p.baseline_evidence,
            clinical_trials=p.clinical_trials,
            similarity_pairs=p.similarity_lookup,
            similarity_threshold=p.similarity_threshold,
            gene_universe=gu,
            con=con,
        )

    results = []
    for t in p.targets_list:
        if return_matched_pairs:
            matched_pairs_df = create_matched_pairs_df(
                genetic_evidence=t,
                clinical_trials=p.clinical_trials,
                similarity_pairs=p.similarity_lookup,
                similarity_threshold=p.similarity_threshold,
                gene_universe=gu,
                con=con,
            )
            target_pairs = set(zip(matched_pairs_df["gene"], matched_pairs_df["ct_efo_id"]))
        else:
            matched_pairs_df = None
            target_pairs = create_matched_pairs_set(
                genetic_evidence=t,
                clinical_trials=p.clinical_trials,
                similarity_pairs=p.similarity_lookup,
                similarity_threshold=p.similarity_threshold,
                gene_universe=gu,
                con=con,
            )

        ct = _prepare_clinical_trials(
            ct=ct_all,
            target_pairs=target_pairs,
            baseline_pairs=baseline_pairs,
            gene_universe=gu,
        )

        enrichment = calculate_all_enrichments(
            df=ct,
            phase_transitions=p.phase_transitions,
            cluster_col=cluster_col,
            bootstrap_replicates=bootstrap_replicates,
            bootstrap_seed=bootstrap_seed,
        )

        results.append(
            _pack_validation_result(
                enrichment,
                ct,
                matched_pairs_df,
                return_trials=return_trials,
                return_matched_pairs=return_matched_pairs,
            ),
        )

    con.close()

    if not batch:
        return results[0]
    return results
