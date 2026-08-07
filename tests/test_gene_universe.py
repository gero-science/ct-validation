"""Tests for gene_universe.build_gene_universe."""

import pandas as pd
import pytest


@pytest.fixture
def target_path(tmp_path):
    """Protein-coding and non-protein-coding targets, including a duplicated symbol
    (two ENSG ids sharing one approvedSymbol) and an Ensembl-ID fallback symbol.
    """
    df = pd.DataFrame(
        {
            "approvedSymbol": [
                "ABCA1",
                "TP53",
                "TP53",
                "MIR21",
                "ENSG00000290734",
            ],
            "biotype": [
                "protein_coding",
                "protein_coding",
                "protein_coding",
                "lncRNA",
                "protein_coding",
            ],
        }
    )
    path = tmp_path / "target.parquet"
    df.to_parquet(path, index=False)
    return path


def test_non_protein_coding_biotypes_are_excluded(gene_universe, target_path, tmp_path):
    """A symbol whose only entry is not protein_coding does not appear in the universe."""
    symbols = gene_universe.build_gene_universe(target_path, tmp_path / "out.txt")

    assert "MIR21" not in symbols


def test_symbol_shared_by_two_ensg_ids_appears_once(gene_universe, target_path, tmp_path):
    """A symbol appearing under two different ENSG ids is written only once."""
    symbols = gene_universe.build_gene_universe(target_path, tmp_path / "out.txt")

    assert symbols.count("TP53") == 1


def test_output_symbols_are_sorted(gene_universe, target_path, tmp_path):
    """Symbols are written in sorted order."""
    symbols = gene_universe.build_gene_universe(target_path, tmp_path / "out.txt")

    assert symbols == sorted(symbols)


def test_ensembl_id_fallback_symbols_are_kept(gene_universe, target_path, tmp_path):
    """A protein-coding gene with no HGNC symbol (approvedSymbol is its Ensembl ID) is
    kept, not filtered out.
    """
    symbols = gene_universe.build_gene_universe(target_path, tmp_path / "out.txt")

    assert "ENSG00000290734" in symbols


def test_output_file_contains_one_symbol_per_line(gene_universe, target_path, tmp_path):
    """The written file matches the returned symbol list, one per line."""
    output_path = tmp_path / "out.txt"
    symbols = gene_universe.build_gene_universe(target_path, output_path)

    assert output_path.read_text().splitlines() == symbols
