#!/usr/bin/env python3
"""Parse ClinVar VCF to gene-disease associations.

Implements King et al. 2019 methodology:
- High-confidence review status filtering
- Pathogenic/Likely pathogenic variants only
- Extracts structured ontology IDs from CLNDISDB (MONDO, Orphanet, MeSH, HP)
- Maps to EFO via direct IDs, OxO mappings, and OLS text search fallback
"""

import argparse
import logging
import time
from collections import defaultdict
from pathlib import Path

import pandas as pd
import requests
import yaml
from cyvcf2 import VCF
from tqdm import tqdm

log = logging.getLogger(__name__)

# High-confidence review status (King et al. 2019)
PATHOGENIC_TERMS = {"Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic"}
HIGH_CONFIDENCE_REVIEW = {
    "criteria_provided,_multiple_submitters,_no_conflicts",
    "reviewed_by_expert_panel",
    "practice_guideline",
}

# Ontology priority for mapping: (column, ontology_name, source_suffix)
ONTOLOGY_PRIORITY = [
    ("mondo_ids", "MONDO", ""),
    ("orphanet_ids", "ORDO", ""),
    ("mesh_ids", "MESH", "_MeSH"),
    ("hp_ids", "HP", ""),
]

# Conditions to exclude from OLS search
UNMAPPABLE_CONDITIONS = {
    "not provided",
    "not_provided",
    "not specified",
    "not_specified",
    "see cases",
    "See_cases",
    "none",
    "inborn_genetic_diseases",
    "inborn genetic diseases",
    "intellectual_disability",
    "intellectual disability",
    "cardiovascular_phenotype",
    "cardiovascular phenotype",
    "hereditary_cancer-predisposing_syndrome",
    "hereditary cancer-predisposing syndrome",
    "neurodevelopmental_disorder",
    "neurodevelopmental disorder",
    "abnormality_of_the_nervous_system",
    "abnormality of the nervous system",
    "abnormality_of_the_musculature",
    "abnormality of the musculature",
}


# === Loading functions ===


def load_efo_ontology(efo_path: Path) -> tuple[set[str], dict[str, str]]:
    """Load EFO ontology IDs and labels from .obo file."""
    efo_ids, efo_labels = set(), {}
    current_id = None

    with efo_path.open() as f:
        for line_ in f:
            line = line_.strip()
            if line.startswith("id:"):
                current_id = line.split()[1].strip()
                efo_ids.add(current_id)
                if current_id.startswith("efo:EFO_"):
                    standard_id = f"EFO:{current_id.replace('efo:EFO_', '')}"
                    efo_ids.add(standard_id)
                    current_id = standard_id
            elif line.startswith("name:") and current_id:
                efo_labels[current_id] = line.replace("name:", "").strip()
                current_id = None

    log.info(f"Loaded {len(efo_ids):,} EFO IDs, {len(efo_labels):,} labels")
    return efo_ids, efo_labels


def load_mapping_file(
    path: Path,
    efo_ids: set[str],
    strip_prefix: str = "",
) -> dict[str, list[str]]:
    """Load ID -> EFO mapping from TSV file."""
    df = pd.read_csv(path, sep="\t")
    df = df[df["mapped_curie"].isin(efo_ids)]

    mapping = defaultdict(list)
    for _, row in df.iterrows():
        source_id = row["curie_id"].replace(strip_prefix, "") if strip_prefix else row["curie_id"]
        if row["mapped_curie"] not in mapping[source_id]:
            mapping[source_id].append(row["mapped_curie"])

    log.info(f"Loaded mapping: {len(mapping):,} IDs from {path.name}")
    return mapping


# === VCF parsing ===


def _extract_clndisdb_ids(clndisdb: str) -> dict[str, list[str]]:
    """Extract ontology IDs from CLNDISDB field."""
    prefixes = {
        "MONDO:": "mondo_ids",
        "Orphanet:": "orphanet_ids",
        "MeSH:": "mesh_ids",
        "HP:": "hp_ids",
    }
    result = {col: [] for col in prefixes.values()}

    for entry in clndisdb.split(","):
        for prefix, col in prefixes.items():
            if prefix in entry:
                result[col].append(prefix + entry.split(prefix)[-1])
                break
    return result


def parse_vcf(vcf_path: Path) -> pd.DataFrame:
    """Parse ClinVar VCF with high-confidence filtering."""
    log.info("Parsing ClinVar VCF...")
    vcf = VCF(str(vcf_path))
    records = []

    for var in tqdm(vcf, desc="Parsing VCF"):
        info = var.INFO

        # Filter by clinical significance
        clnsig = info.get("CLNSIG")
        if not clnsig:
            continue
        sig_terms = {s.strip() for s in str(clnsig).replace("|", ",").split(",") if s.strip()}
        if not (sig_terms & PATHOGENIC_TERMS):
            continue

        # Filter by review status
        clnrevstat = info.get("CLNREVSTAT")
        if not clnrevstat:
            continue
        review_terms = {s.strip() for s in str(clnrevstat).split("|") if s.strip()}
        if not (review_terms & HIGH_CONFIDENCE_REVIEW):
            continue

        # Extract genes and conditions
        clndn = info.get("CLNDN")
        geneinfo = info.get("GENEINFO")
        if not clndn or not geneinfo:
            continue

        genes = [p.split(":")[0] for p in str(geneinfo).split("|") if ":" in p]
        gene_ids = [p.split(":")[1] for p in str(geneinfo).split("|") if ":" in p]
        if not genes:
            continue

        # Extract CLNDISDB ontology IDs
        clndisdb = info.get("CLNDISDB")
        ont_ids = _extract_clndisdb_ids(str(clndisdb)) if clndisdb else {}

        records.append(
            {
                "conditions": str(clndn),
                "gene_symbol": ";".join(genes),
                "gene_id": ";".join(gene_ids),
                **{col: "|".join(ids) if ids else None for col, ids in ont_ids.items()},
            },
        )

    log.info(f"Parsed {len(records):,} high-confidence pathogenic variants")
    return pd.DataFrame(records)


def explode_gene_conditions(df: pd.DataFrame) -> pd.DataFrame:
    """Explode to gene-condition pairs."""
    df["gene"] = df["gene_symbol"].str.split(";")
    df = df.explode("gene")
    df["gene"] = df["gene"].str.strip()

    df["condition"] = df["conditions"].str.split("|")
    df = df.explode("condition")
    df["condition"] = df["condition"].str.strip()
    df = df[df["condition"] != ""]

    log.info(f"Exploded to {len(df):,} gene-condition pairs")
    return df


def apply_condition_filters(df: pd.DataFrame) -> pd.DataFrame:
    """Apply Minikel-style condition filters (somatic, susceptibility, etc.)."""
    exclude_patterns = ["somatic", "susceptibility", "response", r"\?"]
    for pattern in exclude_patterns:
        df = df[~df["condition"].str.contains(pattern, case=False, na=False, regex=True)]
    log.info(f"After condition filters: {len(df):,}")
    return df


# === Condition -> EFO mapping ===


def _find_efo_for_condition(
    row: pd.Series,
    efo_ids: set[str],
    efo_labels: dict[str, str],
    mesh_mapping: dict[str, list[str]],
) -> list[dict] | None:
    """Find EFO mapping for a condition row using priority order."""
    condition = row["condition"]

    for col, ontology, source_suffix in ONTOLOGY_PRIORITY:
        if pd.isna(row.get(col)):
            continue

        for ont_id in row[col].split("|"):
            # MeSH needs mapping lookup, others check direct membership
            if col == "mesh_ids":
                mesh_bare = ont_id.replace("MeSH:", "")
                if mesh_bare in mesh_mapping:
                    return [
                        {
                            "condition": condition,
                            "efo_id": efo_id,
                            "efo_label": efo_labels.get(efo_id, efo_id),
                            "ontology": ontology,
                            "source": f"CLNDISDB{source_suffix}",
                        }
                        for efo_id in mesh_mapping[mesh_bare]
                    ]
            elif ont_id in efo_ids:
                return [
                    {
                        "condition": condition,
                        "efo_id": ont_id,
                        "efo_label": efo_labels.get(ont_id, ont_id),
                        "ontology": ontology,
                        "source": f"CLNDISDB{source_suffix}",
                    },
                ]
    return None


def map_clndisdb(
    df: pd.DataFrame,
    efo_ids: set[str],
    efo_labels: dict[str, str],
    mesh_mapping: dict[str, list[str]],
) -> pd.DataFrame:
    """Map conditions using CLNDISDB IDs (Tier 1)."""
    rows = []
    for _, row in tqdm(df.iterrows(), total=len(df), desc="CLNDISDB mapping"):
        mapping = _find_efo_for_condition(row, efo_ids, efo_labels, mesh_mapping)
        if mapping:
            rows.extend(mapping)

    result = pd.DataFrame(rows)
    if len(result) > 0:
        result = result.drop_duplicates(subset=["condition", "efo_id"], keep="first")
    log.info(f"CLNDISDB mapped {result['condition'].nunique() if len(result) else 0:,} conditions")
    return result


def map_oxo(
    df: pd.DataFrame,
    mapped_conditions: set[str],
    oxo_mapping: dict[str, list[str]],
    efo_labels: dict[str, str],
) -> pd.DataFrame:
    """Map conditions using OxO (Tier 2)."""
    unmapped = df[~df["condition"].isin(mapped_conditions)]
    rows = []

    for _, row in tqdm(unmapped.iterrows(), total=len(unmapped), desc="OxO mapping"):
        condition = row["condition"]
        for col, ontology, _ in ONTOLOGY_PRIORITY:
            if col == "mesh_ids":  # MeSH not in OxO
                continue
            if pd.isna(row.get(col)):
                continue
            for ont_id in row[col].split("|"):
                if ont_id in oxo_mapping:
                    rows.extend(
                        [
                            {
                                "condition": condition,
                                "efo_id": efo_id,
                                "efo_label": efo_labels.get(efo_id, efo_id),
                                "ontology": ontology,
                                "source": f"OxO_{ontology}",
                            }
                            for efo_id in oxo_mapping[ont_id]
                        ],
                    )
                    break
            else:
                continue
            break

    result = pd.DataFrame(rows)
    if len(result) > 0:
        result = result.drop_duplicates(subset=["condition", "efo_id"], keep="first")
    log.info(f"OxO mapped {result['condition'].nunique() if len(result) else 0:,} conditions")
    return result


def map_ols(
    conditions: list[str],
    efo_ids: set[str],
    cache_path: Path | None = None,
) -> pd.DataFrame:
    """Map conditions using OLS text search (Tier 3)."""
    if cache_path and cache_path.exists():
        cached = pd.read_csv(cache_path)
        cached_conditions = set(cached["condition"].unique())
        log.info(f"Loaded {len(cached_conditions):,} cached OLS mappings")
    else:
        cached, cached_conditions = pd.DataFrame(), set()

    to_search = [
        c
        for c in conditions
        if c not in cached_conditions and c.lower() not in UNMAPPABLE_CONDITIONS
    ]
    log.info(f"OLS searching {len(to_search):,} conditions")

    rows = []
    for condition in tqdm(to_search, desc="OLS search"):
        try:
            resp = requests.get(
                "https://www.ebi.ac.uk/ols4/api/search",
                params={
                    "q": condition,
                    "ontology": "efo,mondo",
                    "type": "class",
                    "exact": "false",
                    "rows": 1,
                },
                timeout=10,
            )
            if resp.status_code == 200:
                docs = resp.json().get("response", {}).get("docs", [])
                if docs:
                    doc = docs[0]
                    ont_id = doc.get("obo_id", doc.get("iri", ""))
                    if ont_id and ont_id in efo_ids:
                        rows.append(
                            {
                                "condition": condition,
                                "efo_id": ont_id,
                                "efo_label": doc.get("label", ""),
                                "ontology": doc.get("ontology_name", "").upper(),
                                "source": "OLS_text_search",
                            },
                        )
            time.sleep(0.1)
        except Exception as e:
            log.debug(f"OLS error for '{condition}': {e}")

    result = (
        pd.concat([cached, pd.DataFrame(rows)], ignore_index=True)
        if len(cached)
        else pd.DataFrame(rows)
    )
    if cache_path and len(result) > 0:
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        result.to_csv(cache_path, index=False)

    log.info(f"OLS mapped {result['condition'].nunique() if len(result) else 0:,} conditions")
    return result


# === Main pipeline ===


def parse_clinvar(
    vcf_path: Path,
    efo_path: Path,
    mesh_mapping_path: Path,
    oxo_mapping_path: Path,
    output_path: Path,
    *,
    ols_cache_path: Path | None = None,
    max_efo_per_condition: int = 10,
    skip_ols: bool = False,
) -> pd.DataFrame:
    """Parse ClinVar VCF to gene-disease associations."""
    # Load resources
    efo_ids, efo_labels = load_efo_ontology(efo_path)
    mesh_mapping = load_mapping_file(mesh_mapping_path, efo_ids, strip_prefix="MeSH:")
    oxo_mapping = load_mapping_file(oxo_mapping_path, efo_ids)

    # Parse and filter
    df = parse_vcf(vcf_path)
    df = explode_gene_conditions(df)
    df = apply_condition_filters(df)

    # Tiered mapping
    clndisdb_mapping = map_clndisdb(df, efo_ids, efo_labels, mesh_mapping)
    mapped = set(clndisdb_mapping["condition"].unique()) if len(clndisdb_mapping) else set()

    oxo_result = map_oxo(df, mapped, oxo_mapping, efo_labels)
    mapped |= set(oxo_result["condition"].unique()) if len(oxo_result) else set()

    if skip_ols:
        ols_result = pd.DataFrame()
    else:
        unmapped = df[~df["condition"].isin(mapped)]["condition"].unique().tolist()
        ols_result = map_ols(unmapped, efo_ids, ols_cache_path)
        ols_result = (
            ols_result[~ols_result["condition"].isin(mapped)] if len(ols_result) else ols_result
        )

    final_mapping = pd.concat([clndisdb_mapping, oxo_result, ols_result], ignore_index=True)

    # Filter broad conditions
    efo_counts = final_mapping.groupby("condition")["efo_id"].nunique()
    excessive = efo_counts[efo_counts > max_efo_per_condition].index
    if len(excessive) > 0:
        log.info(
            f"Filtering {len(excessive):,} conditions with >{max_efo_per_condition} EFO mappings",
        )
        final_mapping = final_mapping[~final_mapping["condition"].isin(excessive)]

    # Create gene-disease pairs
    gene_cond = df.groupby(["gene", "condition"], as_index=False).agg(
        gene_id=("gene_id", "first"),
        variant_count=("gene_id", "size"),
    )

    result = (
        gene_cond.merge(final_mapping, on="condition", how="inner")
        .groupby(["gene", "efo_id"], as_index=False)
        .agg(
            {
                "gene_id": "first",
                "variant_count": "sum",
                "efo_label": "first",
                "ontology": "first",
                "source": "first",
            },
        )
        .rename(columns={"variant_count": "n_variants", "efo_label": "trait"})
    )

    log.info(f"Result: {len(result):,} gene-disease pairs, {result['gene'].nunique():,} genes")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--vcf", type=Path, help="ClinVar VCF file")
    parser.add_argument("--efo-ontology", type=Path, help="EFO ontology .obo file")
    parser.add_argument("--mesh-mapping", type=Path, help="MeSH to EFO mapping TSV")
    parser.add_argument("--oxo-mapping", type=Path, help="MONDO/Orphanet/HP to EFO mapping TSV")
    parser.add_argument("--output", type=Path, help="Output parquet path")
    parser.add_argument("--ols-cache", type=Path, help="OLS mapping cache CSV")
    parser.add_argument("--max-efo-per-condition", type=int, help="Max EFO IDs per condition")
    parser.add_argument("--skip-ols", action="store_true", help="Skip OLS text search")
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
        mappings = cfg.get("mappings", {})
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ge_inputs, mappings, thresholds, output_dir = {}, {}, {}, Path()

    parse_clinvar(
        vcf_path=args.vcf or Path(ge_inputs.get("clinvar_vcf", "")),
        efo_path=args.efo_ontology or Path(mappings.get("efo_ontology", "")),
        mesh_mapping_path=args.mesh_mapping or Path(mappings.get("mesh_to_efo", "")),
        oxo_mapping_path=args.oxo_mapping or Path(mappings.get("mondo_orphanet_hp_to_efo", "")),
        output_path=args.output or output_dir / "genetic_evidence" / "clinvar.parquet",
        ols_cache_path=args.ols_cache or output_dir / "genetic_evidence" / "clinvar_ols_cache.csv",
        max_efo_per_condition=args.max_efo_per_condition
        or thresholds.get("clinvar_max_efo_per_condition", 10),
        skip_ols=args.skip_ols,
    )


if __name__ == "__main__":
    main()
