"""Tests for validate() API."""

import warnings

import pandas as pd
import pytest
from ct_validation import MissingOngoingWarning, validate


def test_validate_baseline_mode(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """Baseline mode: YES vs NO groups without baseline exclusion."""
    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    assert isinstance(result, pd.DataFrame)
    assert "rr" in result.columns
    assert len(result) > 0


def test_validate_prioritized_mode(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """Prioritized mode excludes baseline from NO group."""
    baseline_df = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "efo_id": ["EFO:001"],
        }
    )

    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        baseline_evidence=baseline_df,
        similarity_threshold=0.8,
    )

    assert isinstance(result, pd.DataFrame)


def test_validate_args_override_config(tmp_config_yaml, similarity_lookup_df):
    """Explicit args override config values."""
    # Create different data to override config
    custom_ct = pd.DataFrame(
        {
            "gene": ["CUSTOM1"],
            "efo_id": ["EFO:999"],
            "max_phase": [3],
        }
    )
    custom_ge = pd.DataFrame(
        {
            "gene": ["CUSTOM1"],
            "efo_id": ["EFO:999"],
        }
    )
    custom_sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:999"],
            "efo_id_2": ["EFO:999"],
            "similarity": [1.0],
        }
    )

    result = validate(
        config=tmp_config_yaml,
        clinical_trials=custom_ct,
        targets=custom_ge,
        similarity_lookup=custom_sim,
    )

    assert isinstance(result, pd.DataFrame)


def test_validate_returns_dataframe(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """Default return is enrichment DataFrame."""
    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    assert isinstance(result, pd.DataFrame)
    expected_cols = {"phase_from", "phase_to", "phase_label", "x_yes", "n_yes", "rr"}
    assert expected_cols.issubset(set(result.columns))


BOOT_COLS = [
    "rr_boot_ci_lower",
    "rr_boot_ci_upper",
    "or_boot_ci_lower",
    "or_boot_ci_upper",
]


def test_validate_without_cluster_col_has_no_bootstrap_columns(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """Without cluster_col the default output schema is unchanged."""
    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
    )

    assert not set(BOOT_COLS) & set(result.columns)


def test_validate_cluster_col_adds_bootstrap_columns_without_moving_estimates(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """cluster_col adds intervals for both measures and leaves the point estimates alone."""
    kwargs = dict(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
    )
    plain = validate(**kwargs)
    clustered = validate(**kwargs, cluster_col="gene")

    assert set(BOOT_COLS).issubset(set(clustered.columns))
    pd.testing.assert_frame_equal(plain, clustered[plain.columns])


def test_validate_unknown_cluster_col_raises_before_running(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """A column that does not exist is caught up front, not per phase transition."""
    with pytest.raises(KeyError, match="cluster_col"):
        validate(
            clinical_trials=clinical_trials_df,
            targets=genetic_evidence_df,
            similarity_lookup=similarity_lookup_df,
            cluster_col="not_a_column",
        )


def test_validate_bootstrap_seed_is_reproducible(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """Same seed gives identical bounds through validate()."""
    kwargs = dict(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        cluster_col="gene",
        bootstrap_replicates=200,
    )
    first = validate(**kwargs, bootstrap_seed=1)
    repeat = validate(**kwargs, bootstrap_seed=1)

    pd.testing.assert_frame_equal(first, repeat)


def test_validate_return_trials(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """return_trials=True returns tuple of (enrichment, trials)."""
    enrichment, trials = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
        return_trials=True,
    )

    assert isinstance(enrichment, pd.DataFrame)
    assert isinstance(trials, pd.DataFrame)
    assert "has_genetic_evidence" in trials.columns


def test_validate_return_matched_pairs(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """return_matched_pairs=True returns match audit DataFrame."""
    enrichment, matched = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
        return_matched_pairs=True,
    )

    assert isinstance(enrichment, pd.DataFrame)
    assert isinstance(matched, pd.DataFrame)
    assert "p_value" in enrichment.columns
    assert {"gene", "ct_efo_id", "ge_efo_id", "similarity", "match_type"}.issubset(matched.columns)
    assert len(matched) > 0


def test_validate_return_trials_and_matched_pairs(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """Both optional outputs are returned in order."""
    enrichment, trials, matched = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
        return_trials=True,
        return_matched_pairs=True,
    )

    assert isinstance(enrichment, pd.DataFrame)
    assert isinstance(trials, pd.DataFrame)
    assert isinstance(matched, pd.DataFrame)


def test_validate_missing_required_raises():
    """Missing required args raises ValueError."""
    with pytest.raises(ValueError, match="required"):
        validate(clinical_trials=None, targets=None, similarity_lookup=None)


def test_validate_with_gene_universe(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """gene_universe filters to specified genes."""
    gu = {"GENE1", "GENE2"}  # subset of genes

    enrichment, trials = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        gene_universe=gu,
        similarity_threshold=0.8,
        return_trials=True,
    )

    assert trials["gene"].isin(gu).all()


def test_validate_custom_phase_transitions(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """Custom phase_transitions produces correct number of results."""
    custom_transitions = [(1, 2), (1, 4)]

    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
        phase_transitions=custom_transitions,
    )

    assert len(result) == len(custom_transitions)


# --- Batch mode tests ---


def test_validate_batch_returns_list(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """targets=[ge, ge] returns a list of two enrichment DataFrames."""
    results = validate(
        clinical_trials=clinical_trials_df,
        targets=[genetic_evidence_df, genetic_evidence_df],
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    assert isinstance(results, list)
    assert len(results) == 2
    expected_cols = {"phase_from", "phase_to", "phase_label", "x_yes", "n_yes", "rr"}
    for result in results:
        assert isinstance(result, pd.DataFrame)
        assert expected_cols.issubset(set(result.columns))


def test_validate_batch_single_element_list(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """targets=[ge] (one-element list) returns a list of length 1, not unwrapped."""
    results = validate(
        clinical_trials=clinical_trials_df,
        targets=[genetic_evidence_df],
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    assert isinstance(results, list)
    assert len(results) == 1
    assert isinstance(results[0], pd.DataFrame)


def test_validate_batch_return_trials(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """targets=[ge, ge] with return_trials=True returns list of 2 (enrichment, trials) tuples."""
    results = validate(
        clinical_trials=clinical_trials_df,
        targets=[genetic_evidence_df, genetic_evidence_df],
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
        return_trials=True,
    )

    assert isinstance(results, list)
    assert len(results) == 2
    for enrichment, trials in results:
        assert isinstance(enrichment, pd.DataFrame)
        assert isinstance(trials, pd.DataFrame)
        assert "has_genetic_evidence" in trials.columns


def test_validate_batch_matches_single(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """Batch results[0] equals the single-target result for the same inputs."""
    single_result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    batch_results = validate(
        clinical_trials=clinical_trials_df,
        targets=[genetic_evidence_df],
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    pd.testing.assert_frame_equal(batch_results[0], single_result)


# --- similarity_lookup=None (exact matching) tests ---


def test_validate_exact_matching_returns_dataframe(clinical_trials_df, genetic_evidence_df):
    """similarity_lookup=None runs successfully and returns enrichment DataFrame."""
    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=None,
    )

    assert isinstance(result, pd.DataFrame)
    assert len(result) > 0


def test_validate_exact_matching_has_expected_columns(clinical_trials_df, genetic_evidence_df):
    """similarity_lookup=None result contains required enrichment columns."""
    result = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=None,
    )

    expected_cols = {"phase_from", "phase_to", "phase_label", "x_yes", "n_yes", "rr"}
    assert expected_cols.issubset(set(result.columns))


def test_validate_exact_matching_only_counts_exact_pairs(clinical_trials_df):
    """similarity_lookup=None counts only trials with exact (gene, efo_id) GE matches."""
    # GE with only one exact match to the CT data: GENE1/EFO:001
    ge_exact = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})

    # GE with same gene but different efo_id (no CT row exists with EFO:999)
    ge_no_match = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:999"]})

    _, trials_exact = validate(
        clinical_trials=clinical_trials_df,
        targets=ge_exact,
        similarity_lookup=None,
        return_trials=True,
    )
    _, trials_no_match = validate(
        clinical_trials=clinical_trials_df,
        targets=ge_no_match,
        similarity_lookup=None,
        return_trials=True,
    )

    n_matched = trials_exact["has_genetic_evidence"].sum()
    n_unmatched = trials_no_match["has_genetic_evidence"].sum()

    assert n_matched > 0
    assert n_unmatched == 0


def test_validate_exact_matching_fewer_matches_than_similarity(
    clinical_trials_df, genetic_evidence_df, similarity_lookup_df
):
    """Exact matching produces <= matched pairs compared to similarity-expanded matching."""
    _, trials_exact = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=None,
        return_trials=True,
    )
    _, trials_sim = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
        return_trials=True,
    )

    n_exact = trials_exact["has_genetic_evidence"].sum()
    n_sim = trials_sim["has_genetic_evidence"].sum()

    assert n_exact <= n_sim


def test_validate_exact_matching_return_trials(clinical_trials_df, genetic_evidence_df):
    """similarity_lookup=None with return_trials=True returns (enrichment, trials) tuple."""
    enrichment, trials = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=None,
        return_trials=True,
    )

    assert isinstance(enrichment, pd.DataFrame)
    assert isinstance(trials, pd.DataFrame)
    assert "has_genetic_evidence" in trials.columns


def test_validate_exact_matching_gene_universe(clinical_trials_df, genetic_evidence_df):
    """similarity_lookup=None respects gene_universe filtering."""
    gu = {"GENE1", "GENE2"}

    _, trials = validate(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=None,
        gene_universe=gu,
        return_trials=True,
    )

    assert trials["gene"].isin(gu).all()


def test_validate_exact_matching_batch_mode(clinical_trials_df, genetic_evidence_df):
    """similarity_lookup=None works in batch mode and returns a list."""
    results = validate(
        clinical_trials=clinical_trials_df,
        targets=[genetic_evidence_df, genetic_evidence_df],
        similarity_lookup=None,
    )

    assert isinstance(results, list)
    assert len(results) == 2
    for result in results:
        assert isinstance(result, pd.DataFrame)


def test_validate_batch_different_targets(clinical_trials_df, similarity_lookup_df):
    """Batch results differ when target sets have different evidence."""
    ge_large = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5"],
            "efo_id": ["EFO:001", "EFO:010", "EFO:001", "EFO:020", "EFO:002", "EFO:030"],
        }
    )
    ge_small = pd.DataFrame(
        {
            "gene": ["GENE1"],
            "efo_id": ["EFO:001"],
        }
    )

    results = validate(
        clinical_trials=clinical_trials_df,
        targets=[ge_large, ge_small],
        similarity_lookup=similarity_lookup_df,
        similarity_threshold=0.8,
    )

    assert isinstance(results[0], pd.DataFrame)
    assert isinstance(results[1], pd.DataFrame)
    # The two target sets cover different amounts of evidence; rr values must differ
    assert not results[0]["rr"].equals(results[1]["rr"])


@pytest.fixture
def ongoing_clinical_trials_df() -> pd.DataFrame:
    """Two phase-1 pairs: the still-running one has evidence, the concluded one does not,
    so n_yes/n_no say which pair was censored.
    """
    return pd.DataFrame(
        {
            "gene": ["GENE1", "GENE3"],
            "efo_id": ["EFO:001", "EFO:001"],
            "max_phase": [1, 1],
            "is_ongoing": [True, False],
        }
    )


def _phase_row(result: pd.DataFrame, label: str) -> pd.Series:
    return result[result["phase_label"] == label].iloc[0]


def test_validate_censors_ongoing_pairs(ongoing_clinical_trials_df, genetic_evidence_df):
    """is_ongoing on clinical_trials reaches the enrichment counts."""
    result = validate(clinical_trials=ongoing_clinical_trials_df, targets=genetic_evidence_df)

    row = _phase_row(result, "I→II")
    assert (row["n_yes"], row["n_no"]) == (0, 1)  # GENE1 is still running at phase 1


def test_validate_without_ongoing_column_keeps_every_pair(
    ongoing_clinical_trials_df,
    genetic_evidence_df,
):
    """Dropping the column opts out of censoring."""
    result = validate(
        clinical_trials=ongoing_clinical_trials_df.drop(columns=["is_ongoing"]),
        targets=genetic_evidence_df,
    )

    row = _phase_row(result, "I→II")
    assert (row["n_yes"], row["n_no"]) == (1, 1)


@pytest.mark.parametrize(
    "values,dtype",
    [
        ([True, None], "object"),
        ([True, pd.NA], "boolean"),
        ([1.0, float("nan")], "float64"),
    ],
)
def test_missing_ongoing_value_counts_as_concluded(values, dtype, genetic_evidence_df):
    """A missing flag falls back to the no-column default whatever the dtype carrying it."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE3"],
            "efo_id": ["EFO:001", "EFO:001"],
            "max_phase": [1, 1],
            "is_ongoing": pd.Series(values, dtype=dtype),
        }
    )

    result = validate(clinical_trials=ct, targets=genetic_evidence_df)

    row = _phase_row(result, "I→II")
    assert (row["n_yes"], row["n_no"]) == (0, 1)


def test_non_boolean_ongoing_column_raises(genetic_evidence_df):
    """Strings would otherwise blow up with a bare TypeError from inside the counting."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE3"],
            "efo_id": ["EFO:001", "EFO:001"],
            "max_phase": [1, 1],
            "is_ongoing": ["true", "false"],
        }
    )

    with pytest.raises(ValueError, match="must be boolean"):
        validate(clinical_trials=ct, targets=genetic_evidence_df)


def test_validate_rejects_duplicate_gene_indication_pairs(genetic_evidence_df):
    """Repeated pairs would inflate n, so they are refused rather than collapsed."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE1", "GENE2"],
            "efo_id": ["EFO:001", "EFO:001", "EFO:001"],
            "max_phase": [1, 3, 2],
        }
    )

    with pytest.raises(ValueError, match="duplicate"):
        validate(clinical_trials=ct, targets=genetic_evidence_df)


def test_validate_accepts_one_row_per_pair(clinical_trials_df, genetic_evidence_df):
    """The duplicate check does not fire on a properly aggregated input."""
    result = validate(clinical_trials=clinical_trials_df, targets=genetic_evidence_df)

    assert len(result) > 0


def test_validate_warns_when_ongoing_column_is_missing(clinical_trials_df, genetic_evidence_df):
    """A silently uncensored run is announced."""
    with pytest.warns(MissingOngoingWarning, match="is_ongoing"):
        validate(clinical_trials=clinical_trials_df, targets=genetic_evidence_df)


def test_uncensored_run_can_be_silenced_by_category(clinical_trials_df, genetic_evidence_df):
    """Running an uncensored convention on purpose does not have to mean silencing everything."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        warnings.filterwarnings("ignore", category=MissingOngoingWarning)
        warnings.warn("unrelated", UserWarning, stacklevel=1)
        validate(clinical_trials=clinical_trials_df, targets=genetic_evidence_df)

    assert [str(w.message) for w in caught] == ["unrelated"]


def test_validate_warns_once_per_run_not_once_per_target(
    clinical_trials_df,
    genetic_evidence_df,
):
    """The input contract is checked once, not re-checked for every target set."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        validate(
            clinical_trials=clinical_trials_df,
            targets=[genetic_evidence_df, genetic_evidence_df, genetic_evidence_df],
        )

    assert len([w for w in caught if "is_ongoing" in str(w.message)]) == 1


def test_validate_does_not_warn_when_ongoing_column_is_present(
    ongoing_clinical_trials_df,
    genetic_evidence_df,
):
    """No warning when there is nothing to warn about."""
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        validate(clinical_trials=ongoing_clinical_trials_df, targets=genetic_evidence_df)

    assert not [w for w in caught if "is_ongoing" in str(w.message)]
