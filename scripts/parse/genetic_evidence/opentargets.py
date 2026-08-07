#!/usr/bin/env python3
"""Parse OpenTargets genetic evidence datasources to gene-disease associations.

Uses association_by_datasource_direct, restricted to germline variant-disease
datasources; drug/clinical, text-mined (europepmc), model-organism and
functional-genomics sources are excluded to avoid circular evidence with the clinical
trials datasets. uniprot_literature is UniProt-curated variant evidence, not text
mining, so it is kept.

`--indirect` reads association_by_datasource_indirect, the ontology-propagated counterpart
of the same table, and writes to a separate output; every other filter is unchanged, so the
two outputs differ only in propagation.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)

GERMLINE_DATASOURCES = {
    "gwas_credible_sets",
    "eva",
    "gene_burden",
    "genomics_england",
    "gene2phenotype",
    "clingen",
    "uniprot_variants",
    "uniprot_literature",
    "orphanet",
}

# Opt-in: somatic mutations are not germline support for a target. Filtered before the
# groupby, or a somatic score could win the max.
SOMATIC_DATASOURCES = {"eva_somatic", "intogen"}


def load_associations(
    association_dir: Path,
    min_score: float,
    *,
    include_somatic: bool = False,
) -> pd.DataFrame:
    """Load and filter OpenTargets associations to germline genetic datasources."""
    log.info(f"Loading associations from {association_dir}")
    df = pd.read_parquet(
        association_dir,
        columns=["datasourceId", "targetId", "diseaseId", "score"],
    )
    n_total = len(df)

    # Keep only genetic evidence datasources
    datasources = GERMLINE_DATASOURCES | (SOMATIC_DATASOURCES if include_somatic else set())
    df = df[df["datasourceId"].isin(datasources)]
    log.info(
        f"Kept {len(df):,}/{n_total:,} rows from genetic datasources: "
        f"{sorted(df['datasourceId'].unique())}",
    )

    # Aggregate per target-disease pair: max score, collect contributing datasources
    df = df.groupby(["targetId", "diseaseId"], as_index=False).agg(
        score=("score", "max"),
        subsource=("datasourceId", lambda x: "&".join(sorted(set(x)))),
    )

    df = df[df["score"] >= min_score]
    log.info(f"After score filter: {len(df):,} associations (score >= {min_score})")
    return df


def load_target_mapping(target_path: Path) -> dict[str, str]:
    """Load ENSG -> gene symbol mapping from OpenTargets target table."""
    log.info(f"Loading target mapping from {target_path}")
    targets = pd.read_parquet(target_path, columns=["id", "approvedSymbol"])
    mapping = dict(zip(targets["id"], targets["approvedSymbol"]))
    log.info(f"Loaded {len(mapping):,} target mappings")
    return mapping


def parse_opentargets(
    association_dir: Path,
    target_path: Path,
    output_path: Path,
    min_score: float = 0.5,
    *,
    include_somatic: bool = False,
) -> pd.DataFrame:
    """Parse OpenTargets genetic evidence associations to gene-disease pairs."""
    df = load_associations(association_dir, min_score, include_somatic=include_somatic)

    # Map ENSG to gene symbols
    ensg_to_symbol = load_target_mapping(target_path)
    df["gene"] = df["targetId"].map(ensg_to_symbol)
    n_unmapped = df["gene"].isna().sum()
    if n_unmapped:
        log.warning(f"Dropping {n_unmapped:,} rows with unmapped ENSG IDs")
        df = df[df["gene"].notna()]

    # Normalize disease IDs (MONDO_0001234 -> MONDO:0001234)
    df["efo_id"] = df["diseaseId"].str.replace("_", ":")

    # Keep best score per gene-disease pair
    result = df.sort_values("score", ascending=False).drop_duplicates(
        subset=["gene", "efo_id"],
    )

    log.info(f"Result: {len(result):,} gene-disease pairs, {result['gene'].nunique():,} genes")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--association-dir", type=Path, help="OpenTargets association directory")
    parser.add_argument("--target", type=Path, help="OpenTargets target parquet directory")
    parser.add_argument("--output", type=Path, help="Output parquet path")
    parser.add_argument("--min-score", type=float, help="Minimum association score")
    parser.add_argument(
        "--indirect",
        action="store_true",
        help="Read the ontology-propagated association table instead of the direct one",
    )
    parser.add_argument(
        "--include-somatic",
        action="store_true",
        help=f"Also keep somatic datasources ({', '.join(sorted(SOMATIC_DATASOURCES))})",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    if args.indirect and args.association_dir:
        parser.error("--indirect and --association-dir both choose the input; pass only one")
    if args.association_dir and not args.output:
        parser.error("--association-dir requires --output; the default path is the direct table")

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
        thresholds = cfg.get("thresholds", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ge_inputs, thresholds, output_dir = {}, {}, Path()

    association_key = (
        "opentargets_association_indirect" if args.indirect else "opentargets_association"
    )
    output_name = "opentargets_indirect.parquet" if args.indirect else "opentargets.parquet"

    parse_opentargets(
        association_dir=args.association_dir or Path(ge_inputs.get(association_key, "")),
        target_path=args.target or Path(ge_inputs.get("opentargets_target", "")),
        output_path=args.output or output_dir / "genetic_evidence" / output_name,
        min_score=args.min_score
        if args.min_score is not None
        else thresholds.get("opentargets_min_score", 0.5),
        include_somatic=args.include_somatic,
    )


if __name__ == "__main__":
    main()
