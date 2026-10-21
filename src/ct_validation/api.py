"""Main API for clinical trial validation."""

import warnings
from dataclasses import dataclass
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import pyarrow as pa

from ct_validation.config import load_config
from ct_validation.config.schema import DEFAULT_PHASE_TRANSITIONS, Config, Thresholds
from ct_validation.data.schema import CLINICAL_TRIALS, GENETIC_EVIDENCE, SIMILARITY_LOOKUP
from ct_validation.validation.enrichment import calculate_all_enrichments
from ct_validation.validation.matching import (
    _setup_source,
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


class UnmatchableDiseaseWarning(UserWarning):
    """Raised when inputs carry efo_ids that the similarity lookup does not contain.

    Matching joins through the lookup, diagonal included, so such an id can never match:
    clinical-trial pairs carrying one are left out of both groups, and a target row with one
    supports nothing. The causes are a source coded against a different EFO release than the
    lookup (see `efo_remap_obsolete` in configs/parsing.yaml), and terms of ontologies EFO does
    not import, which reach the sources through cross-references. Its own class so callers can
    silence exactly this, rather than every UserWarning their program emits.
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


def _is_arrow_string(dtype) -> bool:
    """
    Tell whether a dtype holds Arrow-backed strings (pandas 3 str, string[pyarrow]).

    Not pd.api.types.is_string_dtype, which also matches categoricals on pandas 2.
    """
    if isinstance(dtype, pd.StringDtype):
        return dtype.storage == "pyarrow"
    if isinstance(dtype, pd.ArrowDtype):
        arrow_type = dtype.pyarrow_dtype
        return pa.types.is_string(arrow_type) or pa.types.is_large_string(arrow_type)
    return False


def _matching_source(data: _DataInput, columns: list[str]) -> _DataInput:
    """
    Narrow a DataFrame to the columns matching reads, with Arrow-backed strings as object.

    duckdb converts every Arrow-backed string column of a registered frame at each query,
    read or not, while it scans object strings in place. An all-null column keeps its
    dtype, so duckdb still types it as a string.
    """
    if not isinstance(data, pd.DataFrame):
        return data
    df = data[columns]
    arrow = [c for c in columns if _is_arrow_string(df[c].dtype) and df[c].notna().any()]
    return df.astype(dict.fromkeys(arrow, object)) if arrow else df


def _lookup_ids(
    similarity_lookup: _DataInput,
    column: str,
    con: duckdb.DuckDBPyConnection,
) -> pd.Series | pd.api.extensions.ExtensionArray:
    """Distinct non-null ids in one column of the unfiltered similarity_lookup."""
    # duckdb would first convert every Arrow string column of the frame; on object
    # strings and paths its DISTINCT is the faster one.
    if isinstance(similarity_lookup, pd.DataFrame) and _is_arrow_string(
        similarity_lookup[column].dtype,
    ):
        return similarity_lookup[column].dropna().unique()
    _setup_source(con, "sim", similarity_lookup)
    return con.execute(
        # sql
        f"SELECT DISTINCT {column} FROM sim WHERE {column} IS NOT NULL",
    ).df()[column]


def _drop_unmatchable_diseases(
    ct: pd.DataFrame,
    similarity_lookup: _DataInput | None,
    con: duckdb.DuckDBPyConnection,
) -> pd.DataFrame:
    """
    Leave out clinical-trial pairs whose efo_id is absent from similarity_lookup.

    Matching joins clinical trials on the lookup's efo_id_2, so such a pair cannot be
    compared with any evidence.
    """
    if similarity_lookup is None:
        return ct

    lookup_ids = _lookup_ids(similarity_lookup, "efo_id_2", con)
    unmatchable = ~ct["efo_id"].isin(lookup_ids)
    if not unmatchable.any():
        return ct

    # vc > 0 skips unused categories of a categorical efo_id.
    vc = ct.loc[unmatchable, "efo_id"].value_counts()
    examples = ", ".join(vc[vc > 0].index[:3])
    if unmatchable.all():
        lookup_examples = ", ".join(sorted(lookup_ids)[:3])
        raise ValueError(
            f"clinical_trials has no efo_id present in similarity_lookup (e.g. {examples}; "
            f"the lookup has {lookup_examples}); check that both use the same id format, "
            "e.g. EFO:0000270 rather than EFO_0000270",
        )

    n_dropped = int(unmatchable.sum())
    warnings.warn(
        f"clinical_trials: {n_dropped:,} of {len(ct):,} pairs ({n_dropped / len(ct):.1%}) "
        f"have an efo_id absent from similarity_lookup (e.g. {examples}) and are left out: "
        "no evidence can match them, so they belong to neither group. Silence this with "
        "warnings.filterwarnings('ignore', category=UnmatchableDiseaseWarning).",
        UnmatchableDiseaseWarning,
        stacklevel=3,
    )
    return ct[~unmatchable]


def _warn_unmatchable(
    con: duckdb.DuckDBPyConnection,
    inputs: dict[str, _DataInput],
    similarity_lookup: _DataInput,
) -> None:
    """Issue an UnmatchableDiseaseWarning per evidence input with efo_ids absent from the lookup."""
    lookup_ids = _lookup_ids(similarity_lookup, "efo_id_1", con)
    con.register("unmatchable_ids", pd.DataFrame({"efo_id": lookup_ids}))
    for i, (name, data) in enumerate(inputs.items()):
        view = f"unmatchable_input_{i}"
        _setup_source(con, view, _matching_source(data, ["efo_id"]))
        # sql
        n, missing, examples = con.execute(
            f"""
            SELECT count(*),
                   count(*) FILTER (WHERE s.efo_id IS NULL),
                   list(DISTINCT d.efo_id ORDER BY d.efo_id) FILTER (WHERE s.efo_id IS NULL)[:3]
            FROM {view} d
            LEFT JOIN unmatchable_ids s ON d.efo_id = s.efo_id
            WHERE d.efo_id IS NOT NULL
            """,
        ).fetchone()
        # Releases the projected copy before matching starts
        con.unregister(view)
        if missing:
            warnings.warn(
                f"{name}: {missing:,} of {n:,} rows ({100 * missing / n:.1f}%) have an efo_id "
                f"absent from similarity_lookup (e.g. {', '.join(examples)}), so they can never "
                "match: terms of another EFO release, or of an ontology EFO does not import. "
                "Silence this with "
                "warnings.filterwarnings('ignore', category=UnmatchableDiseaseWarning).",
                UnmatchableDiseaseWarning,
                stacklevel=3,
            )
    con.unregister("unmatchable_ids")


def _matching_lookup(
    similarity_lookup: _DataInput | None,
    threshold: float,
    con: duckdb.DuckDBPyConnection,
) -> pd.DataFrame | None:
    """
    Load the lookup rows at or above threshold once, for every matching query to read.

    duckdb pushes no filter into a registered frame and plans its joins from the full row
    count. Matching applies the threshold again, so this only narrows what it scans.
    """
    if similarity_lookup is None:
        return None

    columns = ["efo_id_1", "efo_id_2", "similarity"]
    if isinstance(similarity_lookup, pd.DataFrame):
        # A nullable NA becomes NaN and compares False, as NULL does in duckdb.
        similarity = similarity_lookup["similarity"].to_numpy("float64", na_value=np.nan)
        lookup = similarity_lookup.loc[similarity >= threshold]
    else:
        lookup = con.execute(
            # sql
            "SELECT efo_id_1, efo_id_2, similarity FROM read_parquet($path) "
            "WHERE similarity >= $threshold",
            {"path": str(similarity_lookup), "threshold": threshold},
        ).df()
    return _matching_source(lookup, columns)


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

    # .array: pandas 2 fails to build the index from a Series of Arrow dictionary dtype.
    keys = pd.MultiIndex.from_arrays([ct["gene"].array, ct["efo_id"].array])
    has_target = keys.isin(list(target_pairs))
    has_baseline = keys.isin(list(baseline_pairs))

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
            yourself before a validation loop). A clinical_trials pair whose efo_id is
            absent from the lookup cannot be compared with any evidence, so it is left
            out of both groups and an UnmatchableDiseaseWarning is issued. If None,
            performs exact matching on (gene, efo_id) without expansion.
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

    ct_all = _drop_unmatchable_diseases(ct_all, p.similarity_lookup, con)

    if p.similarity_lookup is not None:
        # Clinical trials are covered by _drop_unmatchable_diseases
        inputs: dict[str, _DataInput] = {}
        if p.baseline_evidence is not None:
            inputs["baseline_evidence"] = p.baseline_evidence
        if batch:
            inputs |= {f"targets[{i}]": t for i, t in enumerate(p.targets_list)}
        else:
            inputs["targets"] = p.targets_list[0]
        _warn_unmatchable(con, inputs, p.similarity_lookup)

    similarity_pairs = _matching_lookup(p.similarity_lookup, p.similarity_threshold, con)
    ct_source = _matching_source(p.clinical_trials, ["gene", "efo_id"])

    # Match baseline once (shared across all target sets)
    baseline_pairs: set[tuple[str, str]] = set()
    if p.baseline_evidence is not None:
        baseline_pairs = create_matched_pairs_set(
            genetic_evidence=_matching_source(p.baseline_evidence, ["gene", "efo_id"]),
            clinical_trials=ct_source,
            similarity_pairs=similarity_pairs,
            similarity_threshold=p.similarity_threshold,
            gene_universe=gu,
            con=con,
        )

    # Only the columns enrichment reads, as it copies its input per transition. An
    # absent is_ongoing stays absent: that is what runs it uncensored.
    enrichment_columns = {"max_phase", "has_genetic_evidence", "is_ongoing", cluster_col}

    results = []
    for t in p.targets_list:
        ge_source = _matching_source(t, ["gene", "efo_id"])
        if return_matched_pairs:
            matched_pairs_df = create_matched_pairs_df(
                genetic_evidence=ge_source,
                clinical_trials=ct_source,
                similarity_pairs=similarity_pairs,
                similarity_threshold=p.similarity_threshold,
                gene_universe=gu,
                con=con,
            )
            target_pairs = set(zip(matched_pairs_df["gene"], matched_pairs_df["ct_efo_id"]))
        else:
            matched_pairs_df = None
            target_pairs = create_matched_pairs_set(
                genetic_evidence=ge_source,
                clinical_trials=ct_source,
                similarity_pairs=similarity_pairs,
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
            df=ct.loc[:, ct.columns.isin(enrichment_columns)],
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
