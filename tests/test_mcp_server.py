"""Tests for mcp_server tools and _to_input helper."""

import json

import pandas as pd
import pytest

mcp = pytest.importorskip("mcp", reason="mcp not installed")
from ct_validation.mcp_server import _to_input, ct_validate, expand_disease_set  # noqa: E402

# ---------------------------------------------------------------------------
# _to_input
# ---------------------------------------------------------------------------

SAMPLE_ROWS = [{"gene": "G1", "efo_id": "EFO:001"}, {"gene": "G2", "efo_id": "EFO:002"}]


def test_to_input_file_path_string_returned_as_is():
    path = "/some/path/data.parquet"
    result = _to_input(path, "test")
    assert result == path


def test_to_input_file_path_with_slash_returned_as_is():
    path = "/data/folder/file.csv"
    result = _to_input(path, "test")
    assert result == path


def test_to_input_list_of_dicts_returns_dataframe():
    result = _to_input(SAMPLE_ROWS, "test")
    assert isinstance(result, pd.DataFrame)
    assert list(result.columns) == ["gene", "efo_id"]


def test_to_input_list_of_dicts_has_correct_row_count():
    result = _to_input(SAMPLE_ROWS, "test")
    assert len(result) == 2


def test_to_input_json_string_returns_dataframe():
    json_str = json.dumps(SAMPLE_ROWS)
    result = _to_input(json_str, "test")
    assert isinstance(result, pd.DataFrame)


def test_to_input_json_string_has_correct_columns():
    json_str = json.dumps(SAMPLE_ROWS)
    result = _to_input(json_str, "test")
    assert set(result.columns) == {"gene", "efo_id"}


def test_to_input_none_returns_none():
    assert _to_input(None, "test") is None


def test_to_input_empty_list_raises_value_error():
    with pytest.raises(ValueError, match="empty list"):
        _to_input([], "test")


def test_to_input_empty_json_array_raises_value_error():
    with pytest.raises(ValueError, match="empty list"):
        _to_input("[]", "test")


def test_to_input_label_appears_in_error_message():
    with pytest.raises(ValueError, match="my_label"):
        _to_input([], "my_label")


# ---------------------------------------------------------------------------
# ct_validate — inline data
# ---------------------------------------------------------------------------

CLINICAL_TRIALS_ROWS = [
    {"gene": "GENE1", "efo_id": "EFO:001", "max_phase": 1},
    {"gene": "GENE1", "efo_id": "EFO:002", "max_phase": 2},
    {"gene": "GENE2", "efo_id": "EFO:001", "max_phase": 2},
    {"gene": "GENE2", "efo_id": "EFO:003", "max_phase": 3},
    {"gene": "GENE3", "efo_id": "EFO:001", "max_phase": 1},
    {"gene": "GENE3", "efo_id": "EFO:004", "max_phase": 4},
    {"gene": "GENE4", "efo_id": "EFO:002", "max_phase": 3},
    {"gene": "GENE4", "efo_id": "EFO:005", "max_phase": 2},
    {"gene": "GENE5", "efo_id": "EFO:003", "max_phase": 4},
    {"gene": "GENE5", "efo_id": "EFO:006", "max_phase": 1},
]

TARGETS_ROWS = [
    {"gene": "GENE1", "efo_id": "EFO:001"},
    {"gene": "GENE1", "efo_id": "EFO:010"},
    {"gene": "GENE2", "efo_id": "EFO:001"},
    {"gene": "GENE3", "efo_id": "EFO:020"},
    {"gene": "GENE4", "efo_id": "EFO:002"},
]

SIMILARITY_ROWS = [
    {"efo_id_1": "EFO:001", "efo_id_2": "EFO:001", "similarity": 1.0},
    {"efo_id_1": "EFO:010", "efo_id_2": "EFO:002", "similarity": 0.85},
    {"efo_id_1": "EFO:020", "efo_id_2": "EFO:004", "similarity": 0.75},
    {"efo_id_1": "EFO:001", "efo_id_2": "EFO:005", "similarity": 0.55},
    {"efo_id_1": "EFO:002", "efo_id_2": "EFO:002", "similarity": 1.0},
]


def test_ct_validate_inline_returns_list_of_dicts():
    result = ct_validate(
        clinical_trials=CLINICAL_TRIALS_ROWS,
        targets=TARGETS_ROWS,
        similarity_lookup=SIMILARITY_ROWS,
    )
    assert isinstance(result, list)
    assert len(result) > 0
    assert isinstance(result[0], dict)


def test_ct_validate_inline_has_expected_keys():
    result = ct_validate(
        clinical_trials=CLINICAL_TRIALS_ROWS,
        targets=TARGETS_ROWS,
        similarity_lookup=SIMILARITY_ROWS,
    )
    expected_keys = {
        "phase_from",
        "phase_to",
        "phase_label",
        "n_yes",
        "n_no",
        "x_yes",
        "x_no",
        "rr",
        "rr_ci_lower",
        "rr_ci_upper",
        "p_value",
    }
    assert expected_keys.issubset(result[0].keys())


def test_ct_validate_default_phase_transitions_count():
    result = ct_validate(
        clinical_trials=CLINICAL_TRIALS_ROWS,
        targets=TARGETS_ROWS,
        similarity_lookup=SIMILARITY_ROWS,
    )
    # Default transitions: (1,2), (2,3), (3,4), (1,4)
    assert len(result) == 4


def test_ct_validate_similarity_none_exact_matching():
    result = ct_validate(
        clinical_trials=CLINICAL_TRIALS_ROWS,
        targets=TARGETS_ROWS,
        similarity_lookup=None,
    )
    assert isinstance(result, list)
    assert len(result) > 0


def test_ct_validate_custom_phase_transitions():
    result = ct_validate(
        clinical_trials=CLINICAL_TRIALS_ROWS,
        targets=TARGETS_ROWS,
        similarity_lookup=None,
        phase_transitions=[[1, 2], [1, 4]],
    )
    assert len(result) == 2


def test_ct_validate_custom_phase_transitions_labels():
    result = ct_validate(
        clinical_trials=CLINICAL_TRIALS_ROWS,
        targets=TARGETS_ROWS,
        similarity_lookup=None,
        phase_transitions=[[1, 2], [1, 4]],
    )
    phase_froms = {row["phase_from"] for row in result}
    phase_tos = {row["phase_to"] for row in result}
    assert phase_froms == {1}
    assert phase_tos == {2, 4}


def test_ct_validate_file_paths(tmp_parquet_files):
    result = ct_validate(
        clinical_trials=str(tmp_parquet_files["clinical_trials"]),
        targets=str(tmp_parquet_files["genetic_evidence"]),
        similarity_lookup=str(tmp_parquet_files["similarity"]),
    )
    assert isinstance(result, list)
    assert len(result) > 0


def test_ct_validate_file_paths_has_expected_keys(tmp_parquet_files):
    result = ct_validate(
        clinical_trials=str(tmp_parquet_files["clinical_trials"]),
        targets=str(tmp_parquet_files["genetic_evidence"]),
        similarity_lookup=str(tmp_parquet_files["similarity"]),
    )
    expected_keys = {"phase_from", "phase_to", "rr"}
    assert expected_keys.issubset(result[0].keys())


def test_ct_validate_file_paths_similarity_none(tmp_parquet_files):
    result = ct_validate(
        clinical_trials=str(tmp_parquet_files["clinical_trials"]),
        targets=str(tmp_parquet_files["genetic_evidence"]),
        similarity_lookup=None,
    )
    assert isinstance(result, list)
    assert len(result) > 0


# ---------------------------------------------------------------------------
# expand_disease_set
# ---------------------------------------------------------------------------


def test_expand_disease_set_inline_returns_list():
    result = expand_disease_set(
        efo_ids=["EFO:001"],
        similarity_lookup=SIMILARITY_ROWS,
        similarity_threshold=0.8,
    )
    assert isinstance(result, list)


def test_expand_disease_set_includes_input_id():
    result = expand_disease_set(
        efo_ids=["EFO:001"],
        similarity_lookup=SIMILARITY_ROWS,
        similarity_threshold=0.8,
    )
    # EFO:001 -> EFO:001 with similarity 1.0, so it must appear
    assert "EFO:001" in result


def test_expand_disease_set_includes_similar_disease():
    # EFO:010 -> EFO:002 with similarity 0.85, above threshold 0.8
    result = expand_disease_set(
        efo_ids=["EFO:010"],
        similarity_lookup=SIMILARITY_ROWS,
        similarity_threshold=0.8,
    )
    assert "EFO:002" in result


def test_expand_disease_set_excludes_below_threshold():
    # EFO:020 -> EFO:004 with similarity 0.75, below threshold 0.8
    result = expand_disease_set(
        efo_ids=["EFO:020"],
        similarity_lookup=SIMILARITY_ROWS,
        similarity_threshold=0.8,
    )
    assert "EFO:004" not in result


def test_expand_disease_set_result_is_sorted():
    result = expand_disease_set(
        efo_ids=["EFO:001", "EFO:010"],
        similarity_lookup=SIMILARITY_ROWS,
        similarity_threshold=0.8,
    )
    assert result == sorted(result)


def test_expand_disease_set_file_path(tmp_parquet_files):
    result = expand_disease_set(
        efo_ids=["EFO:001"],
        similarity_lookup=str(tmp_parquet_files["similarity"]),
        similarity_threshold=0.8,
    )
    assert isinstance(result, list)
    assert len(result) > 0


def test_expand_disease_set_threshold_one_returns_only_input():
    # threshold >= 1.0 → no expansion, only input ids (via get_expanded_disease_set logic)
    result = expand_disease_set(
        efo_ids=["EFO:001"],
        similarity_lookup=SIMILARITY_ROWS,
        similarity_threshold=1.0,
    )
    # With threshold 1.0 get_expanded_disease_set returns only the input set
    assert result == ["EFO:001"]
