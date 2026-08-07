"""Tests for the Minikel benchmark-data parser."""

import pandas as pd
import pytest


@pytest.fixture
def write_inputs(tmp_path):
    """Write a pp.tsv and a MeSH-EFO mapping, returning the paths."""

    def _write(programmes: list[dict], mesh_to_efo: list[tuple[str, str]]):
        pp_path = tmp_path / "pp.tsv"
        pd.DataFrame(programmes).to_csv(pp_path, sep="\t", index=False)

        mapping_path = tmp_path / "mesh_to_efo.tsv"
        pd.DataFrame(
            [{"curie_id": f"MeSH:{mesh}", "mapped_curie": efo} for mesh, efo in mesh_to_efo],
        ).to_csv(mapping_path, sep="\t", index=False)

        return pp_path, mapping_path, tmp_path / "out.parquet"

    return _write


def _programme(gene, mesh, ccatnum, acat=None):
    return {
        "ti_uid": f"{gene}-{mesh}",
        "gene": gene,
        "indication_mesh_id": mesh,
        "ccatnum": ccatnum,
        "acat": acat,
    }


def test_citeline_phases_shift_down_by_one(minikel, write_inputs):
    """ccatnum 1=Preclinical..5=Launched maps to max_phase 0=Preclinical..4=Approved."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=5)],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result["max_phase"].tolist() == [4]


def test_active_programme_is_ongoing(minikel, write_inputs):
    """An active category means nothing concluded at the phase it sits in."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=4, acat="Phase III")],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result["is_ongoing"].tolist() == [True]


def test_inactive_programme_is_not_ongoing(minikel, write_inputs):
    """No active category means the programme concluded where it stopped."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=4)],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result["is_ongoing"].tolist() == [False]


def test_concluded_programme_outranks_active_sibling_at_the_same_phase(minikel, write_inputs):
    """A determinate failure at the highest phase beats a sibling still running there."""
    pp, mapping, out = write_inputs(
        [
            _programme("GENE1", "D001", ccatnum=4, acat="Phase III"),
            _programme("GENE1", "D002", ccatnum=4),
        ],
        [("D001", "EFO:001"), ("D002", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert len(result) == 1
    assert result["is_ongoing"].tolist() == [False]


def test_active_programme_above_a_concluded_one_stays_ongoing(minikel, write_inputs):
    """Reaching further keeps the pair undetermined, even with a concluded sibling below."""
    pp, mapping, out = write_inputs(
        [
            _programme("GENE1", "D001", ccatnum=4, acat="Phase III"),
            _programme("GENE1", "D002", ccatnum=3),
        ],
        [("D001", "EFO:001"), ("D002", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result["max_phase"].tolist() == [3]
    assert result["is_ongoing"].tolist() == [True]


def test_one_row_per_gene_indication_pair(minikel, write_inputs):
    """A MeSH term mapping to several EFO terms yields one row per resulting pair."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=3)],
        [("D001", "EFO:001"), ("D001", "MONDO:002")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert sorted(result["efo_id"]) == ["EFO:001", "MONDO:002"]
    assert not result.duplicated(subset=["gene", "efo_id"]).any()


def test_unmapped_indications_are_dropped(minikel, write_inputs):
    """A MeSH term with no EFO cross-reference cannot be matched, so it is left out."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=3), _programme("GENE2", "D999", ccatnum=3)],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result["gene"].tolist() == ["GENE1"]


def test_launched_programmes_are_not_undetermined(minikel, write_inputs):
    """Launched is a terminal category, so an approved pair has a known outcome."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=5, acat="Launched")],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result["max_phase"].tolist() == [4]
    assert result["is_ongoing"].tolist() == [False]


def test_no_mappable_indications_produces_an_empty_result(minikel, write_inputs):
    """An empty result is written rather than raising while computing the summary."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D999", ccatnum=3)],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert result.empty
    assert out.exists()


def test_output_is_accepted_by_validate(minikel, write_inputs):
    """The parquet is directly usable as clinical_trials, censoring included."""
    from ct_validation import validate

    pp, mapping, out = write_inputs(
        [
            _programme("GENE1", "D001", ccatnum=2, acat="Phase I"),
            _programme("GENE2", "D001", ccatnum=2),
        ],
        [("D001", "EFO:001")],
    )
    minikel.parse_clinical_trials(pp, mapping, out)
    targets = pd.DataFrame({"gene": ["GENE1"], "efo_id": ["EFO:001"]})

    result = validate(clinical_trials=out, targets=targets)

    row = result[result["phase_label"] == "I→II"].iloc[0]
    assert (row["n_yes"], row["n_no"]) == (0, 1)  # GENE1 is still running at phase I


def test_indications_without_genetic_insight_are_excluded(minikel, write_inputs, tmp_path):
    """Minikel restrict the denominator to indications where genetic evidence could exist."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=3), _programme("GENE2", "D002", ccatnum=3)],
        [("D001", "EFO:001"), ("D002", "EFO:002")],
    )
    indications = tmp_path / "indic.tsv"
    pd.DataFrame(
        [
            {"indication_mesh_id": "D001", "genetic_insight": "both"},
            {"indication_mesh_id": "D002", "genetic_insight": "none"},
        ],
    ).to_csv(indications, sep="\t", index=False)

    result = minikel.parse_clinical_trials(pp, mapping, out, indications_path=indications)

    assert result["gene"].tolist() == ["GENE1"]


def test_all_indications_kept_without_an_indication_table(minikel, write_inputs):
    """The filter is opt-in, so omitting indic.tsv keeps every indication."""
    pp, mapping, out = write_inputs(
        [_programme("GENE1", "D001", ccatnum=3), _programme("GENE2", "D002", ccatnum=3)],
        [("D001", "EFO:001"), ("D002", "EFO:002")],
    )

    result = minikel.parse_clinical_trials(pp, mapping, out)

    assert sorted(result["gene"]) == ["GENE1", "GENE2"]


@pytest.fixture
def write_associations(tmp_path):
    """Write an assoc.tsv and a MeSH-EFO mapping, returning the paths."""

    def _write(associations: list[dict], mesh_to_efo: list[tuple[str, str]]):
        assoc_path = tmp_path / "assoc.tsv"
        # l2g_share is NaN for the gene-level sources; only OTG rows carry one
        frame = pd.DataFrame(associations)
        if "l2g_share" not in frame.columns:
            frame["l2g_share"] = float("nan")
        frame.to_csv(assoc_path, sep="\t", index=False)

        mapping_path = tmp_path / "mesh_to_efo.tsv"
        pd.DataFrame(
            [{"curie_id": f"MeSH:{mesh}", "mapped_curie": efo} for mesh, efo in mesh_to_efo],
        ).to_csv(mapping_path, sep="\t", index=False)

        return assoc_path, mapping_path, tmp_path / "ge.parquet"

    return _write


def test_associations_map_to_one_row_per_pair(minikel, write_associations):
    """Repeated associations for a pair collapse, and their sources are collected in order."""
    assoc, mapping, out = write_associations(
        [
            {"gene": "GENE1", "mesh_id": "D001", "source": "PICCOLO"},
            {"gene": "GENE1", "mesh_id": "D001", "source": "OTG"},
            {"gene": "GENE1", "mesh_id": "D001", "source": "Genebass"},
        ],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_associations(assoc, mapping, out)

    assert len(result) == 1
    assert result["source"].tolist() == ["Genebass&OTG&PICCOLO"]


def test_unmapped_association_indications_are_dropped(minikel, write_associations):
    """An association whose MeSH term has no EFO cross-reference cannot be matched."""
    assoc, mapping, out = write_associations(
        [
            {"gene": "GENE1", "mesh_id": "D001", "source": "OTG"},
            {"gene": "GENE2", "mesh_id": "D999", "source": "OTG"},
        ],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_associations(assoc, mapping, out)

    assert result["gene"].tolist() == ["GENE1"]


def test_somatic_associations_are_dropped(minikel, write_associations):
    """intOGen is somatic, so it is not germline support for a target."""
    assoc, mapping, out = write_associations(
        [
            {"gene": "GENE1", "mesh_id": "D001", "source": "OTG"},
            {"gene": "GENE2", "mesh_id": "D001", "source": "intOGen"},
        ],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_associations(assoc, mapping, out)

    assert result["gene"].tolist() == ["GENE1"]


def test_locus_to_gene_predictions_below_half_the_locus_score_are_dropped(
    minikel,
    write_associations,
):
    """OTG is a locus-to-gene prediction, so it only supports the gene carrying the locus."""
    assoc, mapping, out = write_associations(
        [
            {"gene": "GENE1", "mesh_id": "D001", "source": "OTG", "l2g_share": 0.5},
            {"gene": "GENE2", "mesh_id": "D001", "source": "OTG", "l2g_share": 0.49},
        ],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_associations(assoc, mapping, out)

    assert result["gene"].tolist() == ["GENE1"]


def test_gene_level_sources_survive_the_locus_score_filter(minikel, write_associations):
    """OMIM and the like are gene-level already, so their missing l2g_share is not a rejection."""
    assoc, mapping, out = write_associations(
        [{"gene": "GENE1", "mesh_id": "D001", "source": "OMIM"}],
        [("D001", "EFO:001")],
    )

    result = minikel.parse_associations(assoc, mapping, out)

    assert result["gene"].tolist() == ["GENE1"]


def test_associations_match_the_targets_schema(minikel, write_associations):
    """The parquet is directly usable as targets."""
    assoc, mapping, out = write_associations(
        [{"gene": "GENE1", "mesh_id": "D001", "source": "OTG"}],
        [("D001", "EFO:001")],
    )

    minikel.parse_associations(assoc, mapping, out)

    assert list(pd.read_parquet(out).columns) == ["gene", "efo_id", "source"]
