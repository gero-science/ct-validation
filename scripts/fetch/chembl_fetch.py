#!/usr/bin/env python3
"""Fetch gene-drug mappings from ChEMBL API.

Queries ChEMBL for human single-protein targets, retrieves activities
above a pChEMBL threshold, and fetches drug info for active compounds.
Saves raw gene-drug mapping as parquet.

This is slow (~hours with 50 parallel jobs). Run once, then use the
output as input to the chembl parser.

Dependencies: chembl_webresource_client, joblib
"""

import argparse
import logging
import time
from pathlib import Path
from typing import Any

import pandas as pd
from chembl_webresource_client.new_client import new_client
from joblib import Parallel, delayed, parallel_config
from tqdm import tqdm

log = logging.getLogger(__name__)

target_chembl = new_client.target
activity_chembl = new_client.activity
drug_chembl = new_client.drug


def get_gene_symbol(target_components: list[dict[str, Any]]) -> str | None:
    """Extract gene symbol from target components."""
    if not target_components:
        return None
    try:
        synonyms = target_components[0].get("target_component_synonyms", [])
        for syn in synonyms:
            if syn.get("syn_type") == "GENE_SYMBOL":
                return syn.get("component_synonym")
    except (IndexError, AttributeError, KeyError):
        pass
    return None


def load_or_fetch_targets(cache_path: Path) -> pd.DataFrame:
    """Load human single protein targets from cache or ChEMBL."""
    if cache_path.exists():
        log.info(f"Loading cached targets from {cache_path}")
        return pd.read_parquet(cache_path)

    log.info("Fetching targets from ChEMBL API...")
    targets = target_chembl.filter(
        target_type="SINGLE PROTEIN",
        organism="Homo sapiens",
    ).only(["organism", "pref_name", "target_chembl_id", "tax_id", "target_components"])

    df = pd.DataFrame(targets)
    df["gene"] = df["target_components"].apply(get_gene_symbol)

    cache_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(cache_path, index=False)
    log.info(f"Cached {len(df):,} targets to {cache_path}")
    return df


def get_activities_for_target(target_chembl_id: str, min_pchembl: float) -> pd.DataFrame | None:
    """Get activities for a target with pChEMBL > threshold."""
    try:
        res = activity_chembl.filter(
            target_chembl_id=target_chembl_id,
            pchembl_value__gt=min_pchembl,
        ).only(["molecule_chembl_id"])
        df = pd.DataFrame(res)
    except Exception as e:
        log.warning(f"Error querying activities for {target_chembl_id}: {e}")
        return None
    else:
        return df if not df.empty else None


def get_drug_info(molecule_ids: list[str]) -> pd.DataFrame | None:
    """Get drug information for molecules."""
    if not molecule_ids:
        return None
    try:
        res = drug_chembl.filter(molecule_chembl_id__in=molecule_ids).only(
            [
                "development_phase",
                "max_phase",
                "first_approval",
                "molecule_chembl_id",
                "synonyms",
                "usan_stem",
                "usan_stem_definition",
                "usan_year",
            ],
        )
        df = pd.DataFrame(res)
        if not df.empty:
            df["max_phase"] = pd.to_numeric(df["max_phase"], errors="coerce")
    except Exception as e:
        log.warning(f"Error querying drug info: {e}")
        return None
    else:
        return df if not df.empty else None


def process_target(target_row: tuple, min_pchembl: float) -> pd.DataFrame | None:
    """Process a single target to get gene-drug mappings."""
    _idx, row = target_row
    target_id = row.target_chembl_id
    gene = row.gene

    if pd.isna(gene):
        return None

    activities = get_activities_for_target(target_id, min_pchembl)
    if activities is None:
        return None

    molecule_ids = activities["molecule_chembl_id"].unique().tolist()
    drugs = get_drug_info(molecule_ids)
    if drugs is None:
        return None

    drugs["initial_target_name"] = gene
    drugs["initial_target_chembl_id"] = target_id

    for col in ["usan_year", "first_approval"]:
        if col in drugs.columns:
            drugs[col] = pd.to_numeric(drugs[col], errors="coerce").astype("Int64")

    return drugs


def fetch_chembl_gene_drug(
    indications_path: Path,
    output: Path,
    min_pchembl: float = 7.0,
    n_jobs: int = 50,
) -> pd.DataFrame:
    """Fetch gene-drug mapping from ChEMBL API."""
    cache_dir = output.parent / "cache"
    cache_dir.mkdir(parents=True, exist_ok=True)

    targets = load_or_fetch_targets(cache_dir / "chembl_targets.parquet")

    log.info(f"Processing {len(targets):,} targets with {n_jobs} jobs...")
    start = time.time()

    with parallel_config(backend="multiprocessing", n_jobs=n_jobs):
        results = Parallel(n_jobs=n_jobs)(
            delayed(process_target)((i, row), min_pchembl)
            for i, row in tqdm(targets.iterrows(), total=len(targets))
        )

    valid = [r for r in results if r is not None]
    log.info(f"Completed in {time.time() - start:.1f}s, {len(valid):,} targets with drugs")

    if not valid:
        log.warning("No gene-drug mappings found")
        return pd.DataFrame()

    df = pd.concat(valid, ignore_index=True)

    # Add molecule names from indications file
    if indications_path.exists():
        indic = pd.read_csv(indications_path, sep="\t")
        names = indic[
            ["Parent Molecule ChEMBL ID", "Parent Molecule Name", "Parent Molecule Type"]
        ].drop_duplicates()
        df = df.merge(
            names,
            left_on="molecule_chembl_id",
            right_on="Parent Molecule ChEMBL ID",
            how="left",
        )
        df = df.drop(columns=["Parent Molecule ChEMBL ID"], errors="ignore")

    output.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(output, index=False)
    log.info(f"Saved {len(df):,} gene-drug records to {output}")
    return df


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--indications", type=Path, required=True, help="ChEMBL indications TSV")
    parser.add_argument("--output", type=Path, required=True, help="Output parquet path")
    parser.add_argument("--min-pchembl", type=float, default=7.0, help="Minimum pChEMBL threshold")
    parser.add_argument("--n-jobs", type=int, default=50, help="Parallel jobs")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    fetch_chembl_gene_drug(args.indications, args.output, args.min_pchembl, args.n_jobs)


if __name__ == "__main__":
    main()
