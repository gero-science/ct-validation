"""Tests for the germline/somatic datasource gate in genetic_evidence/opentargets.py."""

import pandas as pd
import pytest


@pytest.fixture
def associations_path(tmp_path):
    """Three target-disease pairs exercising the somatic gate.

    - TARGET1/DISEASE1: germline (eva) score 0.6 AND somatic (eva_somatic) score 0.9.
      The gate must run before the max-score groupby, so the somatic score must not
      win by default.
    - TARGET2/DISEASE2: somatic-only (intogen) score 0.7. Must not appear by default.
    - TARGET3/DISEASE3: germline (eva) score 0.4, below min_score on its own, AND
      somatic (intogen) score 0.9. Must be dropped entirely by default (not merely
      scored lower), since only the sub-threshold germline evidence counts.
    """
    df = pd.DataFrame(
        {
            "datasourceId": ["eva", "eva_somatic", "intogen", "eva", "intogen"],
            "targetId": ["TARGET1", "TARGET1", "TARGET2", "TARGET3", "TARGET3"],
            "diseaseId": ["DISEASE1", "DISEASE1", "DISEASE2", "DISEASE3", "DISEASE3"],
            "score": [0.6, 0.9, 0.7, 0.4, 0.9],
        }
    )
    path = tmp_path / "associations.parquet"
    df.to_parquet(path, index=False)
    return path


def _pairs(df):
    return set(zip(df["targetId"], df["diseaseId"], strict=True))


def test_somatic_only_pair_absent_by_default(ge_opentargets, associations_path):
    """A pair with no germline evidence at all does not appear by default."""
    result = ge_opentargets.load_associations(associations_path, min_score=0.5)

    assert ("TARGET2", "DISEASE2") not in _pairs(result)


def test_somatic_only_pair_present_with_include_somatic(ge_opentargets, associations_path):
    """The same somatic-only pair appears once somatic sources are opted in."""
    result = ge_opentargets.load_associations(
        associations_path, min_score=0.5, include_somatic=True
    )

    assert ("TARGET2", "DISEASE2") in _pairs(result)


def test_germline_score_used_by_default_even_when_somatic_score_is_higher(
    ge_opentargets, associations_path
):
    """Regression: the somatic gate must run before the max-score groupby.

    TARGET1/DISEASE1 has germline score 0.6 and somatic score 0.9; the somatic score
    must not leak into the max by default.
    """
    result = ge_opentargets.load_associations(associations_path, min_score=0.5)

    row = result[(result["targetId"] == "TARGET1") & (result["diseaseId"] == "DISEASE1")]
    assert row["score"].iloc[0] == pytest.approx(0.6)


def test_somatic_score_wins_when_included(ge_opentargets, associations_path):
    """With somatic sources opted in, the higher somatic score wins the max."""
    result = ge_opentargets.load_associations(
        associations_path, min_score=0.5, include_somatic=True
    )

    row = result[(result["targetId"] == "TARGET1") & (result["diseaseId"] == "DISEASE1")]
    assert row["score"].iloc[0] == pytest.approx(0.9)


def test_pair_dropped_when_only_somatic_evidence_meets_min_score(ge_opentargets, associations_path):
    """A pair whose germline score alone is below min_score is dropped by default,
    even though its excluded somatic score would have cleared the threshold.
    """
    result = ge_opentargets.load_associations(associations_path, min_score=0.5)

    assert ("TARGET3", "DISEASE3") not in _pairs(result)


def test_pair_kept_when_somatic_evidence_included_clears_min_score(
    ge_opentargets, associations_path
):
    """The same pair is kept once somatic evidence is included and clears min_score."""
    result = ge_opentargets.load_associations(
        associations_path, min_score=0.5, include_somatic=True
    )

    assert ("TARGET3", "DISEASE3") in _pairs(result)
