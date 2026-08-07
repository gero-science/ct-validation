"""Shared test fixtures."""

import importlib.util
import sys
from pathlib import Path

import pandas as pd
import pytest

# scripts/parse/clinical_trials/ is a flat script directory (no package, no
# __init__.py); its modules import each other as siblings (e.g. `import
# trial_status` inside opentargets.py). Putting the directory on sys.path once,
# here, lets test modules do the same plain `import trial_status` / `import
# aggregate` and keeps this the single place that knows about that layout.
_CT_SCRIPTS_DIR = Path(__file__).resolve().parent.parent / "scripts" / "parse" / "clinical_trials"
if str(_CT_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_CT_SCRIPTS_DIR))

_PARSE_DIR = Path(__file__).resolve().parent.parent / "scripts" / "parse"
_GENETIC_EVIDENCE_DIR = _PARSE_DIR / "genetic_evidence"


def _load_module_from_path(module_name: str, path: Path):
    """Import a module from an explicit file path, under a private module name.

    scripts/parse/genetic_evidence/ ALSO has its own aggregate.py and opentargets.py
    (distinct from the clinical_trials/ ones already on sys.path above). Adding this
    directory to sys.path too would make `import aggregate` / `import opentargets`
    ambiguous — picking whichever sys.path entry happens to come first — and could
    silently break the clinical-trials tests. Loading by explicit file path sidesteps
    sys.path, and the caller-supplied name keeps these out of sys.modules under a name
    that could collide with the clinical_trials modules.
    """
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def ge_opentargets():
    """The genetic-evidence scripts/parse/genetic_evidence/opentargets.py module."""
    return _load_module_from_path("ge_opentargets", _GENETIC_EVIDENCE_DIR / "opentargets.py")


@pytest.fixture(scope="session")
def gwas_catalog():
    """The scripts/parse/genetic_evidence/gwas_catalog.py module."""
    return _load_module_from_path("ge_gwas_catalog", _GENETIC_EVIDENCE_DIR / "gwas_catalog.py")


@pytest.fixture(scope="session")
def gene_universe():
    """The scripts/parse/gene_universe.py module."""
    return _load_module_from_path("gene_universe", _PARSE_DIR / "gene_universe.py")


@pytest.fixture(scope="session")
def minikel():
    """The scripts/parse/minikel.py module."""
    return _load_module_from_path("minikel", _PARSE_DIR / "minikel.py")


@pytest.fixture
def clinical_trials_df() -> pd.DataFrame:
    """10 clinical trials across phases 1-4."""
    return pd.DataFrame(
        {
            "gene": [
                "GENE1",
                "GENE1",
                "GENE2",
                "GENE2",
                "GENE3",
                "GENE3",
                "GENE4",
                "GENE4",
                "GENE5",
                "GENE5",
            ],
            "efo_id": [
                "EFO:001",
                "EFO:002",
                "EFO:001",
                "EFO:003",
                "EFO:001",
                "EFO:004",
                "EFO:002",
                "EFO:005",
                "EFO:003",
                "EFO:006",
            ],
            "max_phase": [1, 2, 2, 3, 1, 4, 3, 2, 4, 1],
        }
    )


@pytest.fixture
def genetic_evidence_df() -> pd.DataFrame:
    """8 gene-indication associations."""
    return pd.DataFrame(
        {
            "gene": ["GENE1", "GENE1", "GENE2", "GENE3", "GENE4", "GENE5", "GENE6", "GENE7"],
            "efo_id": [
                "EFO:001",
                "EFO:010",
                "EFO:001",
                "EFO:020",
                "EFO:002",
                "EFO:030",
                "EFO:001",
                "EFO:001",
            ],
        }
    )


@pytest.fixture
def similarity_lookup_df() -> pd.DataFrame:
    """EFO similarity pairs (0.5-1.0)."""
    return pd.DataFrame(
        {
            "efo_id_1": ["EFO:001", "EFO:010", "EFO:020", "EFO:030", "EFO:001", "EFO:002"],
            "efo_id_2": ["EFO:001", "EFO:002", "EFO:004", "EFO:003", "EFO:005", "EFO:002"],
            "similarity": [1.0, 0.85, 0.75, 0.6, 0.55, 1.0],
        }
    )


@pytest.fixture
def enrichment_input_df() -> pd.DataFrame:
    """DataFrame ready for enrichment calculation."""
    return pd.DataFrame(
        {
            "max_phase": [1, 1, 2, 2, 3, 3, 4, 4, 1, 2],
            "has_genetic_evidence": [
                True,
                True,
                True,
                False,
                True,
                False,
                True,
                False,
                False,
                False,
            ],
        }
    )


@pytest.fixture
def tmp_parquet_files(tmp_path, clinical_trials_df, genetic_evidence_df, similarity_lookup_df):
    """Write fixtures to parquet files and return paths."""
    ct_path = tmp_path / "clinical_trials.parquet"
    ge_path = tmp_path / "genetic_evidence.parquet"
    sim_path = tmp_path / "similarity.parquet"

    clinical_trials_df.to_parquet(ct_path)
    genetic_evidence_df.to_parquet(ge_path)
    similarity_lookup_df.to_parquet(sim_path)

    return {"clinical_trials": ct_path, "genetic_evidence": ge_path, "similarity": sim_path}


@pytest.fixture
def tmp_config_yaml(tmp_path, tmp_parquet_files):
    """Create a config YAML file pointing to fixture parquets."""
    config_path = tmp_path / "config.yaml"
    config_content = f"""
data:
  clinical_trials: "{tmp_parquet_files["clinical_trials"]}"
  genetic_evidence: "{tmp_parquet_files["genetic_evidence"]}"
  efo_similarity_lookup: "{tmp_parquet_files["similarity"]}"

output:
  dir: "{tmp_path / "output"}"
  save_trials: false

thresholds:
  semantic_similarity: 0.8

phase_transitions:
  - [1, 2]
  - [2, 3]
  - [3, 4]
  - [1, 4]
"""
    config_path.write_text(config_content)
    return config_path
