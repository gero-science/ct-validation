#!/usr/bin/env python3
"""Parse OMIM genemap2 to gene-disease associations.

Implements Minikel-style filtering:
- Mapping code 3 only (established molecular basis)
- Excludes: somatic, susceptibility, response, questionable (?), non-disease ([)
- Cancer susceptibility ({) kept only if explicitly cancer-related
"""

import argparse
import logging
import re
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)

# Minikel filter patterns
EXCLUDE_PATTERNS = {"somatic", "susceptibility", "response"}
CANCER_TERMS = {"cancer", "neoplasm", "tumor", "malignant"}


def parse_phenotypes(phenotype_str: str, gene_symbol: str) -> list[dict]:
    """Parse phenotype string with Minikel-style filtering.

    Format: "Phenotype name, MIMID (mapping_code), inheritance; ..."
    """
    if pd.isna(phenotype_str) or not phenotype_str.strip():
        return []

    results = []
    for pheno_ in phenotype_str.split(";"):
        pheno = pheno_.strip()
        if not pheno:
            continue

        # Extract: phenotype name, MIM ID, mapping code
        match = re.search(r"^(.+?),\s*(\d+)\s*\((\d+)\)", pheno)
        if not match:
            continue

        name, mim_id, mapping_code = match.group(1).strip(), match.group(2), match.group(3)

        # Filter: mapping code 3 only
        if mapping_code != "3":
            continue

        # Filter: questionable (?) or non-disease ([)
        if name.startswith(("?", "[")):
            continue

        name_lower = name.lower()

        # { prefix = susceptibility in OMIM notation; keep only cancer, exclude somatic
        if name.startswith("{"):
            is_cancer = any(t in name_lower for t in CANCER_TERMS)
            if not is_cancer or "somatic" in name_lower:
                continue
        # Filter non-{ entries: somatic, susceptibility, response
        elif any(p in name_lower for p in EXCLUDE_PATTERNS):
            continue

        results.append(
            {
                "gene": gene_symbol,
                "trait": name,
                "phenotype_mim": mim_id,
            },
        )

    return results


def load_genemap(genemap_path: Path) -> pd.DataFrame:
    """Load and parse OMIM genemap2.txt."""
    log.info(f"Loading OMIM genemap2 from {genemap_path}")
    df = pd.read_csv(genemap_path, sep="\t", skiprows=3, dtype=str, na_values=[""])
    df.columns = [col.replace("# ", "") for col in df.columns]

    df = df[df["Phenotypes"].notna() & df["Approved Gene Symbol"].notna()].copy()
    log.info(f"Loaded {len(df):,} entries with phenotypes")
    return df


def load_omim_mapping(mapping_path: Path) -> pd.DataFrame:
    """Load OMIM to EFO/ontology mapping."""
    log.info(f"Loading OMIM mapping from {mapping_path}")
    df = pd.read_csv(mapping_path, sep="\t")
    df["omim_mim"] = df["curie_id"].str.replace("OMIM:", "")
    df = df.rename(columns={"mapped_curie": "ontology_id"})
    log.info(f"Loaded {len(df):,} OMIM mappings")
    return df[["omim_mim", "ontology_id"]]


def parse_omim(
    genemap_path: Path,
    mapping_path: Path,
    output_path: Path,
) -> pd.DataFrame:
    """Parse OMIM to gene-disease associations."""
    # Load data
    genemap = load_genemap(genemap_path)

    # Parse phenotypes
    all_assoc = []
    for _, row in genemap.iterrows():
        all_assoc.extend(parse_phenotypes(row["Phenotypes"], row["Approved Gene Symbol"]))

    df = pd.DataFrame(all_assoc)
    log.info(f"Parsed {len(df):,} gene-phenotype pairs")

    # Map to EFO/ontology IDs
    if mapping_path.exists():
        mapping = load_omim_mapping(mapping_path)
        df = df.merge(mapping, left_on="phenotype_mim", right_on="omim_mim", how="left")
        # Fill unmapped with OMIM: prefix
        df["efo_id"] = df["ontology_id"].fillna("OMIM:" + df["phenotype_mim"])
        n_mapped = df["ontology_id"].notna().sum()
        log.info(f"Mapped {n_mapped:,}/{len(df):,} to ontology IDs")
    else:
        log.warning(f"Mapping file not found: {mapping_path}")
        df["efo_id"] = "OMIM:" + df["phenotype_mim"]

    # Deduplicate
    result = df[["gene", "trait", "efo_id"]].drop_duplicates()
    log.info(
        f"Result: {len(result):,} unique gene-disease pairs, {result['gene'].nunique():,} genes",
    )

    # Save
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--genemap", type=Path, help="OMIM genemap2.txt")
    parser.add_argument("--mapping", type=Path, help="OMIM to EFO mapping TSV")
    parser.add_argument("--output", type=Path, help="Output parquet path")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    config_path = (
        args.config or Path(__file__).parent.parent.parent.parent / "configs" / "parsing.yaml"
    )
    if config_path.exists():
        with config_path.open() as f:
            cfg = yaml.safe_load(f)
        ge_inputs = cfg.get("genetic_evidence", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ge_inputs, output_dir = {}, Path()

    parse_omim(
        genemap_path=args.genemap or Path(ge_inputs.get("omim_genemap2", "")),
        mapping_path=args.mapping or Path(ge_inputs.get("omim_to_efo", "")),
        output_path=args.output or output_dir / "genetic_evidence" / "omim.parquet",
    )


if __name__ == "__main__":
    main()
