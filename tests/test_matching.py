"""Tests for similarity matching."""

import pandas as pd
import pytest
from ct_validation.validation.matching import (
    check_similarity_lookup,
    create_matched_pairs_df,
    create_matched_pairs_set,
    get_expanded_disease_set,
)


def test_exact_match_found():
    """Same gene+efo_id with similarity=1.0 is matched."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:001"], "similarity": [1.0]})

    result = create_matched_pairs_set(ge, ct, sim, similarity_threshold=0.8)

    assert ("GENE1", "EFO:001") in result


def test_similar_match_found():
    """Same gene, similar efo_ids above threshold are matched."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:002"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:002"], "similarity": [0.9]})

    result = create_matched_pairs_set(ge, ct, sim, similarity_threshold=0.8)

    assert ("GENE1", "EFO:002") in result


def test_check_similarity_lookup_rejects_asymmetric():
    """A one-directional (non-symmetric) lookup fails validation loudly."""
    sim = pd.DataFrame({"efo_id_1": ["EFO:002"], "efo_id_2": ["EFO:001"], "similarity": [0.9]})

    with pytest.raises(ValueError, match="not symmetric"):
        check_similarity_lookup(sim, require_diagonal=False)


def test_below_threshold_excluded():
    """Pairs below similarity threshold are not matched."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:002"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:002"], "similarity": [0.7]})

    result = create_matched_pairs_set(ge, ct, sim, similarity_threshold=0.8)

    assert len(result) == 0


def test_different_gene_excluded():
    """Same efo but different gene is not matched."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE2"], "efo_id": ["EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:001"], "similarity": [1.0]})

    result = create_matched_pairs_set(ge, ct, sim, similarity_threshold=0.8)

    assert len(result) == 0


def test_gene_universe_filter():
    """Only genes in universe are matched."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE2"],
            "efo_id": ["EFO:001", "EFO:001"],
            "max_phase": [2, 2],
        }
    )
    ge = pd.DataFrame({"gene": ["GENE1", "GENE2"], "efo_id": ["EFO:001", "EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:001"], "similarity": [1.0]})

    result = create_matched_pairs_set(
        ge, ct, sim, similarity_threshold=0.8, gene_universe={"GENE1"}
    )

    assert ("GENE1", "EFO:001") in result
    assert ("GENE2", "EFO:001") not in result


def test_returns_set_of_tuples():
    """Output is a set of (gene, efo_id) tuples."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:001"], "similarity": [1.0]})

    result = create_matched_pairs_set(ge, ct, sim, similarity_threshold=0.8)

    assert isinstance(result, set)
    assert all(isinstance(p, tuple) and len(p) == 2 for p in result)


def test_empty_inputs_return_empty_set():
    """Empty inputs return empty set without error."""
    ct = pd.DataFrame({"gene": [], "efo_id": [], "max_phase": []})
    ge = pd.DataFrame({"gene": [], "efo_id": []})
    sim = pd.DataFrame({"efo_id_1": [], "efo_id_2": [], "similarity": []})

    result = create_matched_pairs_set(ge, ct, sim, similarity_threshold=0.8)

    assert result == set()


# Tests for get_expanded_disease_set


def test_expand_diseases_exact_match():
    """threshold=1.0 returns only input diseases."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:001"],
            "efo_id_2": ["EFO:001", "EFO:002"],
            "similarity": [1.0, 0.9],
        }
    )

    result = get_expanded_disease_set(["EFO:001"], sim, similarity_threshold=1.0)

    assert result == {"EFO:001"}


def test_expand_diseases_with_similarity():
    """Diseases above threshold are included."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:001"],
            "efo_id_2": ["EFO:002", "EFO:003"],
            "similarity": [0.9, 0.85],
        }
    )

    result = get_expanded_disease_set(["EFO:001"], sim, similarity_threshold=0.8)

    assert result == {"EFO:001", "EFO:002", "EFO:003"}


def test_check_similarity_lookup_accepts_symmetric():
    """A symmetric lookup with a full diagonal passes validation."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:002", "EFO:001", "EFO:002"],
            "efo_id_2": ["EFO:001", "EFO:002", "EFO:002", "EFO:001"],
            "similarity": [1.0, 1.0, 0.9, 0.9],
        }
    )

    # Should not raise.
    assert check_similarity_lookup(sim) is None


def test_expand_diseases_excludes_below_threshold():
    """Diseases below threshold are excluded."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:001"],
            "efo_id_2": ["EFO:002", "EFO:003"],
            "similarity": [0.9, 0.7],
        }
    )

    result = get_expanded_disease_set(["EFO:001"], sim, similarity_threshold=0.8)

    assert "EFO:002" in result
    assert "EFO:003" not in result


def test_expand_diseases_multiple_inputs():
    """Multiple input diseases are all expanded."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:002"],
            "efo_id_2": ["EFO:010", "EFO:020"],
            "similarity": [0.9, 0.9],
        }
    )

    result = get_expanded_disease_set(["EFO:001", "EFO:002"], sim, similarity_threshold=0.8)

    assert result == {"EFO:001", "EFO:002", "EFO:010", "EFO:020"}


def test_expand_diseases_accepts_set():
    """Function accepts set input."""
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:001"], "similarity": [1.0]})

    result = get_expanded_disease_set({"EFO:001"}, sim, similarity_threshold=1.0)

    assert result == {"EFO:001"}


# Tests for similarity_pairs=None (exact matching)


def test_exact_match_none_similarity_found():
    """similarity_pairs=None matches on exact (gene, efo_id)."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})

    result = create_matched_pairs_set(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert ("GENE1", "EFO:001") in result


def test_none_similarity_no_expansion():
    """similarity_pairs=None does not expand to similar EFO terms."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:002"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})

    # EFO:001 and EFO:002 are similar but not identical — no match expected
    result = create_matched_pairs_set(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert len(result) == 0


def test_none_similarity_different_gene_excluded():
    """similarity_pairs=None: same efo_id but different gene is not matched."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE2"], "efo_id": ["EFO:001"]})

    result = create_matched_pairs_set(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert len(result) == 0


def test_none_similarity_multiple_matches():
    """similarity_pairs=None returns all CT pairs with exact GE matches."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE1", "GENE2"],
            "efo_id": ["EFO:001", "EFO:002", "EFO:001"],
            "max_phase": [1, 2, 3],
        }
    )
    ge = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE2", "GENE3"],
            "efo_id": ["EFO:001", "EFO:001", "EFO:099"],
        }
    )

    result = create_matched_pairs_set(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert ("GENE1", "EFO:001") in result
    assert ("GENE2", "EFO:001") in result
    assert ("GENE1", "EFO:002") not in result


def test_none_similarity_returns_set_of_tuples():
    """similarity_pairs=None returns a set of (gene, efo_id) tuples."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})

    result = create_matched_pairs_set(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert isinstance(result, set)
    assert all(isinstance(p, tuple) and len(p) == 2 for p in result)


def test_none_similarity_empty_inputs_return_empty_set():
    """similarity_pairs=None with empty inputs returns empty set."""
    ct = pd.DataFrame({"gene": [], "efo_id": [], "max_phase": []})
    ge = pd.DataFrame({"gene": [], "efo_id": []})

    result = create_matched_pairs_set(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert result == set()


def test_none_similarity_gene_universe_filter():
    """similarity_pairs=None respects gene_universe filtering."""
    ct = pd.DataFrame(
        {
            "gene": ["GENE1", "GENE2"],
            "efo_id": ["EFO:001", "EFO:001"],
            "max_phase": [2, 2],
        }
    )
    ge = pd.DataFrame({"gene": ["GENE1", "GENE2"], "efo_id": ["EFO:001", "EFO:001"]})

    result = create_matched_pairs_set(
        ge, ct, similarity_pairs=None, similarity_threshold=0.8, gene_universe={"GENE1"}
    )

    assert ("GENE1", "EFO:001") in result
    assert ("GENE2", "EFO:001") not in result


# Tests for create_matched_pairs_df


def test_matched_pairs_df_exact_match_columns():
    """Exact matching returns expected audit columns."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})

    result = create_matched_pairs_df(ge, ct, similarity_pairs=None, similarity_threshold=0.8)

    assert list(result.columns) == ["gene", "ct_efo_id", "ge_efo_id", "similarity", "match_type"]
    assert result.iloc[0]["match_type"] == "exact"
    assert result.iloc[0]["similarity"] == pytest.approx(1.0)


def test_matched_pairs_df_similarity_match():
    """Similarity matching records the supporting GE disease ID."""
    ct = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:002"], "max_phase": [2]})
    ge = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})
    sim = pd.DataFrame({"efo_id_1": ["EFO:001"], "efo_id_2": ["EFO:002"], "similarity": [0.9]})

    result = create_matched_pairs_df(ge, ct, sim, similarity_threshold=0.8)

    assert len(result) == 1
    assert result.iloc[0]["ct_efo_id"] == "EFO:002"
    assert result.iloc[0]["ge_efo_id"] == "EFO:001"
    assert result.iloc[0]["match_type"] == "similarity"


def test_check_similarity_lookup_requires_diagonal():
    """Symmetric off-diagonal pairs still fail if the diagonal is missing (default)."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:002"],
            "efo_id_2": ["EFO:002", "EFO:001"],
            "similarity": [0.9, 0.9],
        }
    )

    with pytest.raises(ValueError, match="diagonal"):
        check_similarity_lookup(sim)

    # ...but passes when the diagonal requirement is waived.
    assert check_similarity_lookup(sim, require_diagonal=False) is None


def test_check_similarity_lookup_rejects_weak_diagonal():
    """A diagonal present but below 1.0 is malformed (would drop exact matches)."""
    sim = pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:002", "EFO:001", "EFO:002"],
            "efo_id_2": ["EFO:001", "EFO:002", "EFO:002", "EFO:001"],
            "similarity": [0.5, 1.0, 0.9, 0.9],  # EFO:001 self-similarity < 1.0
        }
    )

    with pytest.raises(ValueError, match="diagonal"):
        check_similarity_lookup(sim)
