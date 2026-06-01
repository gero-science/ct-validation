"""Main API for clinical trial validation."""

from dataclasses import dataclass
from pathlib import Path

import duckdb
import pandas as pd

from ct_validation.config import load_config
from ct_validation.config.schema import DEFAULT_PHASE_TRANSITIONS, Config, Thresholds
from ct_validation.data.schema import CLINICAL_TRIALS, GENETIC_EVIDENCE, SIMILARITY_LOOKUP
from ct_validation.validation.enrichment import calculate_all_enrichments
from ct_validation.validation.matching import create_matched_pairs_df, create_matched_pairs_set

_DataInput = pd.DataFrame | str | Path

# Defaults from config schema (single source of truth)
_DEFAULTS = Thresholds()


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


def _prepare_clinical_trials(
    clinical_trials: pd.DataFrame | str | Path,
    target_pairs: set[tuple[str, str]],
    baseline_pairs: set[tuple[str, str]],
    gene_universe: set[str] | None = None,
) -> pd.DataFrame:
    """
    Load and prepare clinical trials DataFrame for enrichment.

    Returns DataFrame with `has_genetic_evidence` column:
    - True = YES group (has matching targets)
    - False = NO group (no targets AND no baseline)

    Rows with baseline-only support are filtered out.
    """
    if isinstance(clinical_trials, pd.DataFrame):
        ct = clinical_trials.copy()
    else:
        ct = pd.read_parquet(clinical_trials)

    ct = ct[ct["gene"].notna() & ct["efo_id"].notna() & ct["max_phase"].notna()]

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
    return_trials: bool = False,
    return_matched_pairs: bool = False,
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
        clinical_trials: gene, efo_id, max_phase
            Unique gene-indication pairs with max phase reached.
        targets: gene, efo_id
            Gene-indication pairs with supporting genetic evidence.
            Pass a list of DataFrames/paths for batch validation.
        similarity_lookup: efo_id_1, efo_id_2, similarity (optional)
            Symmetric EFO semantic similarity pairs (diagonal has similarity=1.0).
            If None, performs exact matching on (gene, efo_id) without expansion.
        baseline_evidence: gene, efo_id (optional)
            Gene-indication pairs for baseline comparison in prioritized mode.

    Output Schema:
        phase_from, phase_to: Phase transition (e.g., 1→2, 2→3).
        n_yes, n_no: Pairs at phase_from+ (with/without genetic evidence).
        x_yes, x_no: Pairs reaching phase_to+ (with/without genetic evidence).
        rate_yes, rate_no: Progression rates.
        rr, rr_ci_lower, rr_ci_upper: Risk ratio with 95% CI (Katz log method).
        or, or_ci_lower, or_ci_upper: Odds ratio with 95% CI (Woolf logit method).
        p_value: Two-sided Fisher's exact test p-value.

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
    if p.baseline_evidence is not None:
        _validate_input(p.baseline_evidence, GENETIC_EVIDENCE, "baseline_evidence")

    gu = _load_gene_universe(p.gene_universe)

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
            clinical_trials=p.clinical_trials,
            target_pairs=target_pairs,
            baseline_pairs=baseline_pairs,
            gene_universe=gu,
        )

        enrichment = calculate_all_enrichments(
            df=ct,
            phase_transitions=p.phase_transitions,
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
