"""Tests for trial_status.resolve/apply — trial status -> reached/concluded classification."""

import pandas as pd
import pytest
import trial_status


def test_withdrawn_dropped_other_open_targets_statuses_retained():
    """Withdrawn is dropped; every other Open Targets status is retained."""
    status = pd.Series(["Completed", "Withdrawn", "Recruiting", "Withdrawn", "Suspended"])

    keep, _ = trial_status.resolve(
        status, trial_status.OPENTARGETS_CONCLUDED, trial_status.OPENTARGETS_DROPPED
    )

    assert keep.tolist() == [True, False, True, False, True]


def test_withdrawn_dropped_other_trialpanorama_statuses_retained():
    """WITHDRAWN is dropped; every other TrialPanorama status is retained."""
    status = pd.Series(["COMPLETED", "WITHDRAWN", "RECRUITING", "WITHDRAWN"])

    keep, _ = trial_status.resolve(
        status, trial_status.TRIALPANORAMA_CONCLUDED, trial_status.TRIALPANORAMA_DROPPED
    )

    assert keep.tolist() == [True, False, True, False]


@pytest.mark.parametrize("status_value", ["Completed", "Terminated", "Unknown status"])
def test_opentargets_terminal_statuses_are_concluded(status_value):
    """Completed/Terminated/Unknown status are classified as concluded for Open Targets."""
    _, is_concluded = trial_status.resolve(
        pd.Series([status_value]),
        trial_status.OPENTARGETS_CONCLUDED,
        trial_status.OPENTARGETS_DROPPED,
    )

    assert is_concluded.iloc[0]


def test_opentargets_null_status_is_concluded():
    """A NULL status (approval record) is classified as concluded for Open Targets."""
    _, is_concluded = trial_status.resolve(
        pd.Series([None]), trial_status.OPENTARGETS_CONCLUDED, trial_status.OPENTARGETS_DROPPED
    )

    assert is_concluded.iloc[0]


@pytest.mark.parametrize(
    "status_value",
    [
        "Recruiting",
        "Active, not recruiting",
        "Not yet recruiting",
        "Enrolling by invitation",
        "Suspended",
    ],
)
def test_opentargets_open_statuses_are_not_concluded(status_value):
    """Statuses that leave the outcome undetermined are not concluded for Open Targets."""
    _, is_concluded = trial_status.resolve(
        pd.Series([status_value]),
        trial_status.OPENTARGETS_CONCLUDED,
        trial_status.OPENTARGETS_DROPPED,
    )

    assert not is_concluded.iloc[0]


@pytest.mark.parametrize("status_value", ["COMPLETED", "TERMINATED", "UNKNOWN"])
def test_trialpanorama_terminal_statuses_are_concluded(status_value):
    """COMPLETED/TERMINATED/UNKNOWN are classified as concluded for TrialPanorama."""
    _, is_concluded = trial_status.resolve(
        pd.Series([status_value]),
        trial_status.TRIALPANORAMA_CONCLUDED,
        trial_status.TRIALPANORAMA_DROPPED,
    )

    assert is_concluded.iloc[0]


@pytest.mark.parametrize(
    "status_value",
    [
        "RECRUITING",
        "ACTIVE_NOT_RECRUITING",
        "NOT_YET_RECRUITING",
        "ENROLLING_BY_INVITATION",
        "SUSPENDED",
    ],
)
def test_trialpanorama_open_statuses_are_not_concluded(status_value):
    """Statuses that leave the outcome undetermined are not concluded for TrialPanorama."""
    _, is_concluded = trial_status.resolve(
        pd.Series([status_value]),
        trial_status.TRIALPANORAMA_CONCLUDED,
        trial_status.TRIALPANORAMA_DROPPED,
    )

    assert not is_concluded.iloc[0]


def test_unclassified_status_raises():
    """A status in neither the concluded nor the dropped mapping raises, not silently 'ongoing'."""
    status = pd.Series(["Completed", "Some New ClinicalTrials.gov Status"])

    with pytest.raises(ValueError, match="unclassified trial status"):
        trial_status.resolve(
            status, trial_status.OPENTARGETS_CONCLUDED, trial_status.OPENTARGETS_DROPPED
        )


def test_trialpanorama_null_status_raises():
    """TrialPanorama has no approval records, so a NULL status is unclassified and raises."""
    status = pd.Series(["COMPLETED", None])

    with pytest.raises(ValueError, match="unclassified trial status"):
        trial_status.resolve(
            status, trial_status.TRIALPANORAMA_CONCLUDED, trial_status.TRIALPANORAMA_DROPPED
        )


def test_apply_drops_withdrawn_records():
    """apply() removes Withdrawn records and keeps every other row."""
    df = pd.DataFrame({"status": ["Completed", "Withdrawn", "Recruiting"], "id": [1, 2, 3]})

    result = trial_status.apply(
        df, trial_status.OPENTARGETS_CONCLUDED, trial_status.OPENTARGETS_DROPPED
    )

    assert result["id"].tolist() == [1, 3]


def test_apply_is_concluded_is_plain_bool_dtype_with_no_missing_values():
    """apply() adds a real bool column (not nullable/object), never NA."""
    df = pd.DataFrame({"status": ["Completed", "Recruiting", "Terminated"]})

    result = trial_status.apply(
        df, trial_status.OPENTARGETS_CONCLUDED, trial_status.OPENTARGETS_DROPPED
    )

    assert result["is_concluded"].dtype == bool
    assert result["is_concluded"].notna().all()


def test_apply_is_concluded_stays_aligned_with_its_row_after_filtering():
    """Each retained row's is_concluded value corresponds to that row, not a shifted one."""
    df = pd.DataFrame({"status": ["Recruiting", "Withdrawn", "Completed"], "id": ["a", "b", "c"]})

    result = trial_status.apply(
        df, trial_status.OPENTARGETS_CONCLUDED, trial_status.OPENTARGETS_DROPPED
    )

    assert dict(zip(result["id"], result["is_concluded"], strict=True)) == {"a": False, "c": True}
