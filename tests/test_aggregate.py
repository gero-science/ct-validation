"""Tests for aggregate.aggregate_clinical_trials — the end-to-end aggregation pipeline."""

import aggregate
import pandas as pd
import pytest

OUTPUT_TABLES = [
    "drug_id_mapping",
    "gene_drug_mapping",
    "drug_indication_mapping",
    "gene_indication_max_phase",
    "gene_drug_by_source",
    "drug_indication_by_source",
]


def _read_outputs(output_dir) -> dict[str, pd.DataFrame]:
    return {name: pd.read_parquet(output_dir / f"{name}.parquet") for name in OUTPUT_TABLES}


@pytest.fixture(scope="module")
def phase_logic_source_dir(tmp_path_factory):
    """Three sources exercising the is_ongoing/is_concluded phase logic and source exclusion.

    opentargets (carries is_concluded):
      - CHEMBLA/EFO:ONGOING reaches phase 2 with only an ongoing record.
      - CHEMBLA/EFO:TIE reaches phase 3 with both a concluded AND an ongoing
        record at that same phase (concluded must win).
      - CHEMBLB/EFO:CONCLUDED reaches phase 4 with only a concluded record.
      - CHEMBLB/EFO:APPROVED_ONGOING reaches phase 4 with only an ongoing record
        (approval is terminal, so it still counts as determined).
      - CHEMBLB also has a record with no efo_id at all.
    dgidb (stands in for an included source with no is_concluded column; the real
    DGIdb parser emits no indication file):
      - DrugC/EFO:NOCOL reaches phase 3.
    chembl (excluded from the drug-indication leg — it carries no trial status):
      - CHEMBLD/EFO:EXCLUDED reaches phase 4.
      - DrugD-alias/CHEMBL0-ALT appears in no other file and joins DrugD's cluster
        only via the shared drugbank_id DBD; it sorts below CHEMBLD, so the cluster's
        representative chembl_id shows whether it merged. Realistic because only
        ChEMBL's gene-drug file is filtered to human targets with pChEMBL, so its
        indications carry drugs that file lacks.
    """
    source_dir = tmp_path_factory.mktemp("phase_logic_source")

    opentargets_gene_drug = pd.DataFrame(
        {
            "gene": ["GENEA", "GENEB"],
            "drug_name": ["DrugA", "DrugB"],
            "chembl_id": ["CHEMBLA", "CHEMBLB"],
        }
    )
    opentargets_drug_indication = pd.DataFrame(
        {
            "drug_name": ["DrugA", "DrugA", "DrugA", "DrugB", "DrugB", "DrugB"],
            "chembl_id": ["CHEMBLA", "CHEMBLA", "CHEMBLA", "CHEMBLB", "CHEMBLB", "CHEMBLB"],
            "efo_id": [
                "EFO:ONGOING",
                "EFO:TIE",
                "EFO:TIE",
                "EFO:CONCLUDED",
                "EFO:APPROVED_ONGOING",
                None,
            ],
            "phase": [2, 3, 3, 4, 4, 4],
            "is_concluded": [False, True, False, True, False, True],
        }
    )
    dgidb_gene_drug = pd.DataFrame({"gene": ["GENEC"], "drug_name": ["DrugC"]})
    dgidb_drug_indication = pd.DataFrame(
        {
            "drug_name": ["DrugC"],
            "efo_id": ["EFO:NOCOL"],
            "phase": [3],
        }
    )
    chembl_gene_drug = pd.DataFrame(
        {
            "gene": ["GENED"],
            "drug_name": ["DrugD"],
            "chembl_id": ["CHEMBLD"],
            "drugbank_id": ["DBD"],
        }
    )
    chembl_drug_indication = pd.DataFrame(
        {
            "drug_name": ["DrugD", "DrugD-alias"],
            "chembl_id": ["CHEMBLD", "CHEMBL0-ALT"],
            "drugbank_id": ["DBD", "DBD"],
            "efo_id": ["EFO:EXCLUDED", "EFO:EXCLUDED"],
            "phase": [4, 4],
        }
    )

    opentargets_gene_drug.to_parquet(source_dir / "opentargets_gene_drug.parquet", index=False)
    opentargets_drug_indication.to_parquet(
        source_dir / "opentargets_drug_indication.parquet", index=False
    )
    dgidb_gene_drug.to_parquet(source_dir / "dgidb_gene_drug.parquet", index=False)
    dgidb_drug_indication.to_parquet(source_dir / "dgidb_drug_indication.parquet", index=False)
    chembl_gene_drug.to_parquet(source_dir / "chembl_gene_drug.parquet", index=False)
    chembl_drug_indication.to_parquet(source_dir / "chembl_drug_indication.parquet", index=False)

    return source_dir


@pytest.fixture(scope="module")
def phase_logic_result(phase_logic_source_dir, tmp_path_factory):
    """Aggregate phase_logic_source_dir once and return its output tables."""
    output_dir = tmp_path_factory.mktemp("phase_logic_output")
    aggregate.aggregate_clinical_trials(phase_logic_source_dir, output_dir)
    return _read_outputs(output_dir)


@pytest.fixture
def multi_id_cluster_source_dir(tmp_path):
    """One drug reachable under two different chembl_ids via a shared drugbank_id.

    Union-Find must merge CHEMBL1 and CHEMBL2 into a single cluster (they share
    drugbank_id DB1), and drug-indication rows keyed on either chembl_id must
    both survive into the aggregated output.
    """
    source_dir = tmp_path / "source"
    source_dir.mkdir()

    gene_drug = pd.DataFrame(
        {
            "gene": ["GENEX", "GENEX"],
            "drug_name": ["DrugX-A", "DrugX-B"],
            "chembl_id": ["CHEMBL1", "CHEMBL2"],
            "drugbank_id": ["DB1", "DB1"],
        }
    )
    drug_indication = pd.DataFrame(
        {
            "drug_name": ["DrugX-A", "DrugX-B"],
            "chembl_id": ["CHEMBL1", "CHEMBL2"],
            "efo_id": ["EFO:FROM1", "EFO:FROM2"],
            "phase": [2, 3],
        }
    )

    gene_drug.to_parquet(source_dir / "opentargets_gene_drug.parquet", index=False)
    drug_indication.to_parquet(source_dir / "opentargets_drug_indication.parquet", index=False)

    return source_dir


@pytest.fixture(scope="module")
def two_source_result(tmp_path_factory):
    """One drug both sources describe, disagreeing on phase and on action type.

    Open Targets and TrialPanorama cluster together on the normalized drug name.
    """
    source_dir = tmp_path_factory.mktemp("two_source_source")

    pd.DataFrame(
        {
            "gene": ["GENEE"],
            "drug_name": ["DrugE"],
            "chembl_id": ["CHEMBLE"],
            "action_type": ["INHIBITOR"],
        }
    ).to_parquet(source_dir / "opentargets_gene_drug.parquet", index=False)
    pd.DataFrame(
        {
            "drug_name": ["DrugE"],
            "chembl_id": ["CHEMBLE"],
            "efo_id": ["EFO:SHARED"],
            "phase": [2],
            "is_concluded": [True],
        }
    ).to_parquet(source_dir / "opentargets_drug_indication.parquet", index=False)

    pd.DataFrame(
        {"gene": ["GENEE"], "drug_name": ["DrugE"], "action_type": ["AGONIST"]}
    ).to_parquet(source_dir / "trialpanorama_gene_drug.parquet", index=False)
    # Still running at 4, so phase_concluded has to come from Open Targets' 2 while phase
    # comes from here — the two columns resolving to different sources on one pair.
    pd.DataFrame(
        {"drug_name": ["DrugE"], "efo_id": ["EFO:SHARED"], "phase": [4], "is_concluded": [False]}
    ).to_parquet(source_dir / "trialpanorama_drug_indication.parquet", index=False)

    output_dir = tmp_path_factory.mktemp("two_source_output")
    aggregate.aggregate_clinical_trials(source_dir, output_dir)
    return _read_outputs(output_dir)


# --- Per-source tables ---


def test_by_source_keeps_each_source_own_phase(two_source_result):
    """The combined table reports 4; the by-source table still knows OT said 2."""
    assert two_source_result["drug_indication_mapping"]["phase"].tolist() == [4]

    by_source = two_source_result["drug_indication_by_source"]
    assert dict(zip(by_source["source"], by_source["phase"], strict=True)) == {
        "OpenTargets": 2,
        "TrialPanorama": 4,
    }


def test_by_source_keeps_each_source_own_action_type(two_source_result):
    """Action types are `&`-joined when collapsed, so neither can be attributed."""
    assert two_source_result["gene_drug_mapping"]["action_type"].tolist() == ["AGONIST&INHIBITOR"]

    by_source = two_source_result["gene_drug_by_source"]
    assert dict(zip(by_source["source"], by_source["action_type"], strict=True)) == {
        "OpenTargets": "INHIBITOR",
        "TrialPanorama": "AGONIST",
    }


@pytest.mark.parametrize("results", ["phase_logic_result", "two_source_result"])
def test_by_source_tables_neither_invent_nor_lose_edges(results, request):
    """Rolling the per-source rows back up reproduces the combined tables' sources.

    Both fixtures, because in phase_logic_result every pair comes from a single source, so
    on its own it never exercises a roll-up across sources.
    """
    phase_logic_result = request.getfixturevalue(results)
    for combined, by_source, keys in [
        ("gene_drug_mapping", "gene_drug_by_source", ["gene", "our_drug_id"]),
        ("drug_indication_mapping", "drug_indication_by_source", ["our_drug_id", "efo_id"]),
    ]:
        rolled = (
            phase_logic_result[by_source]
            .groupby(keys)["source"]
            .agg(lambda sources: "&".join(sorted(set(sources))))
        )
        expected = phase_logic_result[combined].set_index(keys)["source"].sort_index()
        pd.testing.assert_series_equal(rolled, expected)


@pytest.mark.parametrize("results", ["phase_logic_result", "two_source_result"])
@pytest.mark.parametrize("column", ["phase", "phase_concluded"])
def test_drug_indication_by_source_phases_max_to_the_combined_ones(results, column, request):
    """Both phase columns are exactly the max over the per-source values.

    phase_concluded matters as much as phase: the two together decide is_ongoing, so
    anything reconstructing censoring from the per-source table depends on both. Only
    two_source_result puts two sources on one pair, which is the case worth testing —
    there Open Targets reports phase 2 concluded and TrialPanorama phase 4 still running,
    so the combined row must be phase 4 with phase_concluded 2.

    Dtypes are not compared: a source contributing only ongoing records has a null
    phase_concluded, so the per-source column is nullable where the combined one, holding
    the max across sources, need not be.
    """
    result = request.getfixturevalue(results)
    rolled = result["drug_indication_by_source"].groupby(["our_drug_id", "efo_id"])[column].max()
    expected = (
        result["drug_indication_mapping"].set_index(["our_drug_id", "efo_id"])[column].sort_index()
    )
    pd.testing.assert_series_equal(rolled, expected, check_dtype=False)


# --- Determinism ---


def test_aggregation_is_deterministic_across_runs(phase_logic_source_dir, tmp_path):
    """Running aggregation twice on identical inputs produces identical output.

    Row order is not part of the contract (nothing ORDER BYs the final tables),
    so each output is sorted on its key columns before comparing.
    """
    out1, out2 = tmp_path / "out1", tmp_path / "out2"
    aggregate.aggregate_clinical_trials(phase_logic_source_dir, out1)
    aggregate.aggregate_clinical_trials(phase_logic_source_dir, out2)

    result1, result2 = _read_outputs(out1), _read_outputs(out2)
    for name in OUTPUT_TABLES:
        df1, df2 = result1[name], result2[name]
        key_cols = [c for c in ("our_drug_id", "gene", "efo_id", "source") if c in df1.columns]
        df1 = df1.sort_values(key_cols).reset_index(drop=True)
        df2 = df2.sort_values(key_cols).reset_index(drop=True)
        pd.testing.assert_frame_equal(df1, df2)


# --- Multi-id cluster: silent data loss regression ---


def test_drug_indication_rows_survive_under_every_clustered_id(
    multi_id_cluster_source_dir, tmp_path
):
    """A drug-indication row keyed on ANY id in a merged cluster must not be dropped.

    Before the fix, only the arbitrarily-elected representative identifier's
    rows joined against the drug ID mapping, and the rest of the cluster's rows
    silently vanished from the output.
    """
    output_dir = tmp_path / "output"
    aggregate.aggregate_clinical_trials(multi_id_cluster_source_dir, output_dir)

    drug_indication = pd.read_parquet(output_dir / "drug_indication_mapping.parquet")

    assert drug_indication["our_drug_id"].nunique() == 1
    assert set(drug_indication["efo_id"]) == {"EFO:FROM1", "EFO:FROM2"}


# --- is_ongoing / is_concluded phase logic ---


def test_is_ongoing_matches_its_definition_on_every_row(phase_logic_result):
    """is_ongoing == below approval AND (max_phase_concluded is null OR < max_phase)."""
    gim = phase_logic_result["gene_indication_max_phase"]

    undetermined = gim["max_phase_concluded"].isna() | (
        gim["max_phase_concluded"] < gim["max_phase"]
    )
    assert (gim["is_ongoing"] == ((gim["max_phase"] < 4) & undetermined)).all()


def test_max_phase_concluded_never_exceeds_max_phase(phase_logic_result):
    """max_phase_concluded <= max_phase always, where max_phase_concluded is set."""
    gim = phase_logic_result["gene_indication_max_phase"]

    concluded = gim.dropna(subset=["max_phase_concluded"])
    assert (concluded["max_phase_concluded"] <= concluded["max_phase"]).all()


def test_max_phase_is_never_null(phase_logic_result):
    """max_phase is always populated."""
    gim = phase_logic_result["gene_indication_max_phase"]

    assert gim["max_phase"].notna().all()


def test_pair_with_only_ongoing_records_at_max_phase_is_flagged_ongoing(phase_logic_result):
    """A pair whose highest phase has only ongoing records is flagged is_ongoing."""
    gim = phase_logic_result["gene_indication_max_phase"]

    row = gim.loc[(gim["gene"] == "GENEA") & (gim["efo_id"] == "EFO:ONGOING")].iloc[0]
    assert row["is_ongoing"]


def test_pair_with_concluded_record_at_max_phase_is_not_ongoing(phase_logic_result):
    """A pair with a concluded record at its highest phase is not flagged ongoing."""
    gim = phase_logic_result["gene_indication_max_phase"]

    row = gim.loc[(gim["gene"] == "GENEB") & (gim["efo_id"] == "EFO:CONCLUDED")].iloc[0]
    assert not row["is_ongoing"]


def test_approved_pair_is_never_ongoing(phase_logic_result):
    """Approval is terminal, so reaching phase 4 settles the outcome."""
    gim = phase_logic_result["gene_indication_max_phase"]

    row = gim.loc[(gim["gene"] == "GENEB") & (gim["efo_id"] == "EFO:APPROVED_ONGOING")].iloc[0]
    assert row["max_phase"] == 4
    assert not row["is_ongoing"]


def test_concluded_record_wins_over_ongoing_record_at_same_phase(phase_logic_result):
    """A concluded record at the highest phase overrides a co-occurring ongoing one."""
    gim = phase_logic_result["gene_indication_max_phase"]

    row = gim.loc[(gim["gene"] == "GENEA") & (gim["efo_id"] == "EFO:TIE")].iloc[0]
    assert not row["is_ongoing"]


def test_source_without_is_concluded_column_is_treated_as_concluded(phase_logic_result):
    """A source lacking is_concluded entirely is treated as concluded, not censored."""
    gim = phase_logic_result["gene_indication_max_phase"]

    row = gim.loc[(gim["gene"] == "GENEC") & (gim["efo_id"] == "EFO:NOCOL")].iloc[0]
    assert not row["is_ongoing"]


@pytest.fixture
def unknown_status_source_dir(tmp_path):
    """One pair whose only record at its top phase carries Unknown status."""
    source_dir = tmp_path / "unknown_source"
    source_dir.mkdir()
    pd.DataFrame(
        {"gene": ["GENEU"], "drug_name": ["DrugU"], "chembl_id": ["CHEMBLU"]},
    ).to_parquet(source_dir / "opentargets_gene_drug.parquet", index=False)
    pd.DataFrame(
        {
            "drug_name": ["DrugU"],
            "chembl_id": ["CHEMBLU"],
            "efo_id": ["EFO:UNKNOWN"],
            "phase": [2],
            "status": ["Unknown status"],
            # what the parser emits for Unknown under the default reading
            "is_concluded": [True],
        },
    ).to_parquet(source_dir / "opentargets_drug_indication.parquet", index=False)
    return source_dir


@pytest.mark.parametrize(("unknown_ongoing", "expected"), [(False, False), (True, True)])
def test_unknown_status_sensitivity_arm(
    unknown_status_source_dir,
    tmp_path,
    unknown_ongoing,
    expected,
):
    """Unknown reads as concluded by default; the sensitivity arm censors it instead."""
    output_dir = tmp_path / f"out_{unknown_ongoing}"
    aggregate.aggregate_clinical_trials(
        unknown_status_source_dir,
        output_dir,
        unknown_ongoing=unknown_ongoing,
    )

    gim = pd.read_parquet(output_dir / "gene_indication_max_phase.parquet")
    assert gim["is_ongoing"].tolist() == [expected]


# --- Source exclusion and unmapped indications ---


def test_chembl_indications_are_excluded_from_the_evidence_leg(phase_logic_result):
    """ChEMBL carries no trial status, so its indications must not reach the output.

    Its drugs still contribute identity, so GENED survives in gene_drug_mapping.
    """
    gim = phase_logic_result["gene_indication_max_phase"]

    assert "EFO:EXCLUDED" not in set(gim["efo_id"])
    assert "GENED" not in set(gim["gene"])
    assert "GENED" in set(phase_logic_result["gene_drug_mapping"]["gene"])


def test_chembl_indications_still_contribute_drug_identity(phase_logic_result):
    """Excluding ChEMBL evidence must not drop its drugs from the clustering.

    CHEMBL0-ALT reaches DrugD's cluster only if chembl_drug_indication is still fed
    to _build_drug_id_mapping.
    """
    drug_id_mapping = phase_logic_result["drug_id_mapping"]

    row = drug_id_mapping.loc[drug_id_mapping["drug_name"] == "DrugD"].iloc[0]
    assert row["chembl_id"] == "CHEMBL0-ALT", "CHEMBL0-ALT did not join DrugD's cluster"


def test_indications_without_an_efo_id_produce_no_pair(phase_logic_result):
    """A drug-indication record with no efo_id is not a pair and must be dropped.

    Grouping them yields one row per gene whose max_phase spans every unmapped
    indication it has, inflating any published count.
    """
    gim = phase_logic_result["gene_indication_max_phase"]

    assert gim["efo_id"].notna().all()
    assert "EFO:CONCLUDED" in set(gim.loc[gim["gene"] == "GENEB", "efo_id"])


# --- Join-induced row multiplication ---


def test_gene_indication_max_phase_is_unique_per_gene_efo_pair(phase_logic_result):
    """(gene, efo_id) is unique in gene_indication_max_phase."""
    gim = phase_logic_result["gene_indication_max_phase"]

    assert not gim.duplicated(subset=["gene", "efo_id"]).any()


def test_drug_indication_mapping_is_unique_per_drug_efo_pair(phase_logic_result):
    """(our_drug_id, efo_id) is unique in drug_indication_mapping."""
    drug_indication = phase_logic_result["drug_indication_mapping"]

    assert not drug_indication.duplicated(subset=["our_drug_id", "efo_id"]).any()
