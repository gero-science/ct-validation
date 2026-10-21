"""Tests for validate() API."""

import warnings

import numpy as np
import pandas as pd
import pyarrow as pa
import pytest
from ct_validation import MissingOngoingWarning, UnmatchableDiseaseWarning, validate


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


# --- UnmatchableDiseaseWarning: CT rows whose efo_id is absent from the lookup ---

_LOOKUP_KW = {"similarity_threshold": 0.8}


def _unmatchable(caught):
    return [w for w in caught if issubclass(w.category, UnmatchableDiseaseWarning)]


def _validate_recording(**kwargs):
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        result = validate(**kwargs)
    return result, _unmatchable(caught)


def test_unmatchable_row_is_left_out_like_a_hand_removed_row(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """A CT row with an id absent from the lookup counts as if it were not in the input."""
    kept = clinical_trials_df[clinical_trials_df["efo_id"] != "EFO:006"]

    full, _ = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )
    by_hand, _ = _validate_recording(
        clinical_trials=kept,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    pd.testing.assert_frame_equal(full, by_hand)


def test_unmatchable_warning_names_count_and_example(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """One warning reports how many rows were dropped and an example id."""
    _, caught = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    assert len(caught) == 1
    assert "1 of 10" in str(caught[0].message)
    assert "EFO:006" in str(caught[0].message)


def test_no_unmatchable_warning_when_all_ids_are_in_lookup(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """Nothing to drop, nothing to warn about."""
    ct = clinical_trials_df[clinical_trials_df["efo_id"] != "EFO:006"]

    _, caught = _validate_recording(
        clinical_trials=ct,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    assert caught == []


def test_no_lookup_means_no_drop_and_no_unmatchable_warning(
    clinical_trials_df,
    genetic_evidence_df,
):
    """Exact matching keeps every CT row, including ids no lookup would contain."""
    (_, trials), caught = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=None,
        return_trials=True,
    )

    assert caught == []
    assert len(trials) == len(clinical_trials_df)


def test_unmatchable_warning_once_per_batch_call(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """Batch mode checks the CT once, not once per target set."""
    _, caught = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=[genetic_evidence_df, genetic_evidence_df, genetic_evidence_df],
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    assert len(caught) == 1


def test_returned_trials_exclude_unmatchable_row(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """return_trials shows the rows that were actually counted."""
    (_, trials), _ = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        return_trials=True,
        **_LOOKUP_KW,
    )

    assert "EFO:006" not in set(trials["efo_id"])
    assert len(trials) == len(clinical_trials_df) - 1


def test_unmatchable_drop_works_with_lookup_parquet_path(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
    tmp_parquet_files,
):
    """A lookup given as a parquet path behaves like the DataFrame."""
    from_df, _ = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )
    from_path, caught = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=tmp_parquet_files["similarity"],
        **_LOOKUP_KW,
    )

    pd.testing.assert_frame_equal(from_df, from_path)
    assert len(caught) == 1


def test_every_ct_id_absent_from_lookup_raises(genetic_evidence_df, similarity_lookup_df):
    """A wholesale id-format mismatch raises rather than returning an empty result."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE2"],
            "efo_id": ["EFO_001", "EFO_002"],
            "max_phase": [1, 2],
        }
    )

    with pytest.raises(ValueError, match="similarity_lookup"):
        validate(
            clinical_trials=ct,
            targets=genetic_evidence_df,
            similarity_lookup=similarity_lookup_df,
            **_LOOKUP_KW,
        )


def test_unmatchable_warning_points_at_caller(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """stacklevel attributes the warning to the code that called validate()."""
    _, caught = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    assert caught[0].filename == __file__


def test_unmatchable_targets_are_announced_per_input(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """A target efo_id missing from the lookup supports nothing, so each input says how many."""
    ct = clinical_trials_df[clinical_trials_df["efo_id"] != "EFO:006"]
    ge = pd.concat([genetic_evidence_df, pd.DataFrame({"gene": ["GENE9"], "efo_id": ["EFO:OLD"]})])

    _, caught = _validate_recording(
        clinical_trials=ct,
        targets=[genetic_evidence_df, ge],
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    assert [str(w.message)[:24] for w in caught] == ["targets[1]: 1 of 9 rows "]
    assert "EFO:OLD" in str(caught[0].message)


def test_target_id_only_in_lookup_efo_id_2_is_unmatchable(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """Evidence joins on efo_id_1, so an id the lookup has only as efo_id_2 cannot match."""
    ct = clinical_trials_df[clinical_trials_df["efo_id"] != "EFO:006"]
    ge = pd.concat([genetic_evidence_df, pd.DataFrame({"gene": ["GENE9"], "efo_id": ["EFO:004"]})])

    _, caught = _validate_recording(
        clinical_trials=ct,
        targets=ge,
        similarity_lookup=similarity_lookup_df,
        **_LOOKUP_KW,
    )

    assert len(caught) == 1
    assert "EFO:004" in str(caught[0].message)


# --- Same results for every input representation (Arrow strings, categoricals, paths) ---

_ID_COLS = ("gene", "efo_id", "efo_id_1", "efo_id_2")
_BASELINE = pd.DataFrame({"gene": ["GENE3", "GENE4"], "efo_id": ["EFO:001", "EFO:005"]})
_PAIR_KW = {
    "phase_transitions": [(1, 2), (2, 3)],
    "return_trials": True,
    "return_matched_pairs": True,
}

_STRING_DTYPES = [
    pytest.param("string[pyarrow]", id="string[pyarrow]"),
    pytest.param(pd.ArrowDtype(pa.string()), id="ArrowDtype-string"),
    pytest.param(
        "str",
        id="pandas3-str",
        marks=pytest.mark.skipif(
            int(pd.__version__.split(".")[0]) < 3,
            reason="default str dtype needs pandas >= 3",
        ),
    ),
]


def _with_dtype(df, dtype):
    """Cast the id columns of df, whichever of them it has, to dtype."""
    return df.astype({c: dtype for c in _ID_COLS if c in df.columns})


def _run_quiet(**kwargs):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return validate(**kwargs)


def _objectify(df):
    """Id columns as object: result frames keep whatever string dtype the input had."""
    ids = ("gene", "efo_id", "ct_efo_id", "ge_efo_id")
    return df.astype({c: object for c in ids if c in df.columns})


def _assert_same_batch(got, expected, pairs_as_object=False):
    """Enrichment and matched pairs exactly equal; trials equal once ids are object."""
    assert len(got) == len(expected)
    for (g_enrich, g_trials, g_pairs), (e_enrich, e_trials, e_pairs) in zip(got, expected):
        pd.testing.assert_frame_equal(g_enrich, e_enrich, check_exact=True)
        if pairs_as_object:  # categorical ids pass through to the pairs frame
            g_pairs, e_pairs = _objectify(g_pairs), _objectify(e_pairs)
        pd.testing.assert_frame_equal(g_pairs, e_pairs, check_exact=True)
        pd.testing.assert_frame_equal(_objectify(g_trials), _objectify(e_trials), check_exact=True)


def _batch_inputs(clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """Object-string inputs, the reference: on pandas 3 the fixtures default to str."""
    ge = _with_dtype(genetic_evidence_df, object)
    return {
        "clinical_trials": _with_dtype(clinical_trials_df, object),
        "targets": [ge, ge.iloc[:3]],
        "similarity_lookup": _with_dtype(similarity_lookup_df, object),
        "baseline_evidence": _with_dtype(_BASELINE, object),
        **_PAIR_KW,
    }


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_arrow_string_inputs_give_same_result_as_object_inputs(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    inputs = _batch_inputs(clinical_trials_df, genetic_evidence_df, similarity_lookup_df)
    arrow_inputs = {
        **inputs,
        "clinical_trials": _with_dtype(inputs["clinical_trials"], dtype),
        "targets": [_with_dtype(t, dtype) for t in inputs["targets"]],
        "similarity_lookup": _with_dtype(inputs["similarity_lookup"], dtype),
        "baseline_evidence": _with_dtype(_BASELINE, dtype),
    }

    _assert_same_batch(_run_quiet(**arrow_inputs), _run_quiet(**inputs))


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_missing_values_in_arrow_targets_and_lookup_match_object_result(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    ge = pd.concat(
        [genetic_evidence_df, pd.DataFrame({"gene": [None], "efo_id": ["EFO:001"]})],
        ignore_index=True,
    )
    sim = pd.concat(
        [
            similarity_lookup_df,
            pd.DataFrame(
                {
                    "efo_id_1": [None, "EFO:001"],
                    "efo_id_2": ["EFO:001", None],
                    "similarity": [0.9, 0.9],
                },
            ),
        ],
        ignore_index=True,
    )
    inputs = {
        **_batch_inputs(clinical_trials_df, ge, sim),
        "targets": [_with_dtype(ge, object)],
    }
    arrow_inputs = {
        **inputs,
        "targets": [_with_dtype(ge, dtype)],
        "similarity_lookup": _with_dtype(sim, dtype),
    }

    _assert_same_batch(_run_quiet(**arrow_inputs), _run_quiet(**inputs))


def _unmatchable_message(**kwargs):
    with pytest.warns(UnmatchableDiseaseWarning) as caught:
        validate(**kwargs)
    return str(_unmatchable(caught)[0].message)


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_unmatchable_warning_text_is_the_same_for_arrow_and_object_lookup(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    kwargs = {
        "clinical_trials": clinical_trials_df,
        "targets": genetic_evidence_df,
        **_LOOKUP_KW,
    }

    arrow = _unmatchable_message(
        **kwargs,
        similarity_lookup=_with_dtype(similarity_lookup_df, dtype),
    )
    obj = _unmatchable_message(
        **kwargs,
        similarity_lookup=_with_dtype(similarity_lookup_df, object),
    )

    assert arrow == obj


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_all_unmatchable_error_text_is_the_same_for_arrow_and_object_lookup(
    dtype,
    genetic_evidence_df,
    similarity_lookup_df,
):
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE2", "GENE3"],
            "efo_id": ["EFO_003", "EFO_001", "EFO_002"],
            "max_phase": [1, 2, 3],
        },
    )
    lookup = pd.concat(
        [
            similarity_lookup_df,
            pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": [None], "similarity": [0.9]}),
        ],
        ignore_index=True,
    )
    kwargs = {"clinical_trials": ct, "targets": genetic_evidence_df, **_LOOKUP_KW}

    with pytest.raises(ValueError, match="similarity_lookup") as arrow:
        validate(**kwargs, similarity_lookup=_with_dtype(lookup, dtype))
    with pytest.raises(ValueError, match="similarity_lookup") as obj:
        validate(**kwargs, similarity_lookup=_with_dtype(lookup, object))

    assert str(arrow.value) == str(obj.value)


def _unmatchable_messages(clinical_trials, evidence_inputs):
    _, caught = _validate_recording(
        clinical_trials=clinical_trials,
        **evidence_inputs,
        **_LOOKUP_KW,
    )
    return [str(w.message) for w in caught]


def _warning_inputs(genetic_evidence_df, similarity_lookup_df, dtype):
    """Targets with an absent id and an efo_id_2-only id, a baseline with an absent id."""
    ge = _with_dtype(genetic_evidence_df, object)
    old = pd.DataFrame({"gene": ["GENE9"], "efo_id": ["EFO:OLD"]})
    only_2 = pd.DataFrame({"gene": ["GENE9"], "efo_id": ["EFO:004"]})
    baseline = pd.concat([_BASELINE, old], ignore_index=True)
    return {
        "targets": [
            _with_dtype(pd.concat([ge, t], ignore_index=True), dtype) for t in (old, only_2)
        ],
        "baseline_evidence": _with_dtype(baseline, dtype),
        "similarity_lookup": _with_dtype(similarity_lookup_df, dtype),
    }


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_unmatchable_warnings_are_the_same_for_arrow_and_object_inputs(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    ct = _with_dtype(clinical_trials_df, object)

    arrow = _unmatchable_messages(
        ct,
        _warning_inputs(genetic_evidence_df, similarity_lookup_df, dtype),
    )
    obj = _unmatchable_messages(
        ct,
        _warning_inputs(genetic_evidence_df, similarity_lookup_df, object),
    )

    assert len(obj) == 4  # clinical_trials, baseline_evidence, targets[0], targets[1]
    assert arrow == obj


def test_unmatchable_warnings_are_the_same_for_parquet_and_dataframe_lookup(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
    tmp_path,
):
    path = tmp_path / "similarity.parquet"
    similarity_lookup_df.to_parquet(path)
    inputs = _warning_inputs(genetic_evidence_df, similarity_lookup_df, object)

    from_df = _unmatchable_messages(clinical_trials_df, inputs)
    from_path = _unmatchable_messages(clinical_trials_df, {**inputs, "similarity_lookup": path})

    assert from_path == from_df


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_missing_target_id_is_not_counted_in_arrow_unmatchable_warning(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    ct = _with_dtype(clinical_trials_df, object)
    extra = pd.DataFrame({"gene": ["GENE9", "GENE8"], "efo_id": ["EFO:OLD", None]})
    ge = pd.concat([_with_dtype(genetic_evidence_df, object), extra], ignore_index=True)

    def messages(d):
        return _unmatchable_messages(
            ct,
            {
                "targets": _with_dtype(ge, d),
                "similarity_lookup": _with_dtype(similarity_lookup_df, d),
            },
        )

    obj = messages(object)
    assert "targets: 1 of 9 rows" in obj[-1]
    assert messages(dtype) == obj


def test_float32_similarity_at_threshold_gives_same_result_for_dataframe_and_parquet(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
    tmp_path,
):
    """The GENE1/EFO:002 pair hinges on a float32(0.9) row at threshold 0.9."""
    lookup = similarity_lookup_df.assign(
        similarity=np.array([1.0, 0.9, 0.75, 0.6, 0.55, 1.0], dtype="float32"),
    )
    path = tmp_path / "similarity.parquet"
    lookup.to_parquet(path)
    kwargs = {
        "clinical_trials": clinical_trials_df,
        "targets": genetic_evidence_df,
        "similarity_threshold": 0.9,
        **_PAIR_KW,
    }

    from_df = _run_quiet(**kwargs, similarity_lookup=lookup)
    from_path = _run_quiet(**kwargs, similarity_lookup=path)

    for got, expected in zip(from_df, from_path):
        pd.testing.assert_frame_equal(got, expected, check_exact=True)


@pytest.mark.parametrize("dtype", _STRING_DTYPES)
def test_empty_arrow_target_in_batch_gives_string_typed_matched_pair_ids(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    results = _run_quiet(
        clinical_trials=_with_dtype(clinical_trials_df, dtype),
        targets=[
            _with_dtype(genetic_evidence_df, dtype),
            _with_dtype(genetic_evidence_df.iloc[:0], dtype),
        ],
        similarity_lookup=_with_dtype(similarity_lookup_df, dtype),
        return_matched_pairs=True,
    )

    empty_pairs = results[1][1]
    assert empty_pairs.empty
    assert pd.api.types.is_string_dtype(empty_pairs["ge_efo_id"])


_CATEGORICAL_DTYPES = [
    pytest.param("category", id="categorical"),
    pytest.param(pd.ArrowDtype(pa.dictionary(pa.int8(), pa.string())), id="arrow-dictionary"),
]


@pytest.mark.parametrize("dtype", _CATEGORICAL_DTYPES)
def test_categorical_trial_ids_give_same_result_as_object(
    dtype,
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    inputs = _batch_inputs(clinical_trials_df, genetic_evidence_df, similarity_lookup_df)
    categorical_inputs = {
        **inputs,
        "clinical_trials": clinical_trials_df.astype({"gene": dtype, "efo_id": dtype}),
    }

    _assert_same_batch(
        _run_quiet(**categorical_inputs),
        _run_quiet(**inputs),
        pairs_as_object=True,
    )


def test_trial_id_present_in_lookup_only_below_threshold_stays_in_no_group(
    clinical_trials_df,
    genetic_evidence_df,
    similarity_lookup_df,
):
    """EFO:003 occurs in the lookup only at 0.6 (threshold 0.8): kept, never matched."""
    (enrichment, trials), caught = _validate_recording(
        clinical_trials=clinical_trials_df,
        targets=genetic_evidence_df,
        similarity_lookup=similarity_lookup_df,
        return_trials=True,
        **_LOOKUP_KW,
    )

    below = trials[trials["efo_id"] == "EFO:003"]
    assert len(below) == 2
    assert not below["has_genetic_evidence"].any()
    assert "EFO:003" not in str(caught[0].message)


# --- A parquet input after a DataFrame input is matched on its own rows ---

_FOUR_PAIRS_CT = pd.DataFrame(
    {
        "gene": ["A", "B", "C", "D"],
        "efo_id": ["EFO:1", "EFO:1", "EFO:1", "EFO:1"],
        "max_phase": [4, 1, 4, 1],
    },
)
_FOUR_PAIRS_KW = {"clinical_trials": _FOUR_PAIRS_CT, "phase_transitions": [(1, 4)]}
_ONE_ID_LOOKUP = pd.DataFrame({"efo_id_1": ["EFO:1"], "efo_id_2": ["EFO:1"], "similarity": [1.0]})


def _pairs(*genes):
    return pd.DataFrame({"gene": list(genes), "efo_id": ["EFO:1"] * len(genes)})


@pytest.mark.parametrize("similarity_lookup", [None, _ONE_ID_LOOKUP])
def test_path_target_after_dataframe_target_matches_like_path_target_alone(
    tmp_path,
    similarity_lookup,
):
    path_target = tmp_path / "b.parquet"
    _pairs("B", "C").to_parquet(path_target)
    kw = {**_FOUR_PAIRS_KW, "similarity_lookup": similarity_lookup}

    mixed = _run_quiet(targets=[_pairs("A"), path_target], **kw)
    alone = _run_quiet(targets=path_target, **kw)

    pd.testing.assert_frame_equal(mixed[1], alone, check_exact=True)


def test_path_target_with_dataframe_baseline_matches_dataframe_target(tmp_path):
    path_target = tmp_path / "t.parquet"
    target = _pairs("B", "C")
    target.to_parquet(path_target)
    kw = {**_FOUR_PAIRS_KW, "baseline_evidence": _pairs("A", "D")}

    via_path = _run_quiet(targets=path_target, **kw)
    via_df = _run_quiet(targets=target, **kw)

    pd.testing.assert_frame_equal(via_path, via_df, check_exact=True)
