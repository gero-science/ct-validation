#!/usr/bin/env python3
"""Parse Minikel et al. (2024) benchmark data to the validate() schemas.

pp.tsv (Pharmaprojects programmes) -> clinical trials; assoc.tsv.gz (genetic associations)
-> targets. Both legs are MeSH-coded upstream and are remapped to EFO here.

Citeline phases (ccatnum) are 1=Preclinical..5=Launched; ours are 0=Preclinical..4=Approved.
A programme is still running when Pharmaprojects gives it an active category (acat), which
is how Minikel censor: their per-transition succ_* flags are NA exactly for active
programmes that have not reached the target phase.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)

# Somatic drivers are not germline support for a target; genetic_evidence/opentargets.py
# gates intogen the same way, and Minikel's own published association set excludes it.
SOMATIC_SOURCES = {"intOGen"}

# Open Targets Genetics associations are locus-to-gene predictions, so they only support the
# gene carrying at least half the locus score. Other sources are gene-level, with no l2g_share.
MIN_L2G_SHARE = 0.5


def _optional(path: str | None) -> Path | None:
    """Config paths are optional; an absent key means the step is skipped."""
    return Path(path) if path else None


def load_mesh_to_efo(mapping_path: Path) -> pd.DataFrame:
    """Load MeSH -> EFO cross-references, keyed by bare MeSH ID."""
    log.info(f"Loading MeSH to EFO mappings from {mapping_path}")
    mappings = pd.read_csv(mapping_path, sep="\t")
    mappings["mesh_id"] = mappings["curie_id"].str.removeprefix("MeSH:")
    mappings = mappings[["mesh_id", "mapped_curie"]].drop_duplicates()
    log.info(f"Loaded {len(mappings):,} MeSH-EFO edges")
    return mappings


def parse_clinical_trials(
    pharmaprojects_path: Path,
    mesh_to_efo_path: Path,
    output_path: Path,
    indications_path: Path | None = None,
) -> pd.DataFrame:
    """Parse Pharmaprojects programmes to censored gene-indication pairs."""
    log.info(f"Loading Pharmaprojects from {pharmaprojects_path}")
    pp = pd.read_csv(pharmaprojects_path, sep="\t", low_memory=False)
    log.info(f"Loaded {len(pp):,} programmes")

    pp = pp[pp["gene"].notna() & pp["indication_mesh_id"].notna() & pp["ccatnum"].notna()]

    if indications_path is not None:
        # Minikel restrict the denominator to indications where genetic evidence could exist
        # at all, and their published RS is computed after it: their reported cohort sizes
        # match this filter to within 8 programmes and differ from the unfiltered ones by
        # up to 1,748. Dropping it would make this arm no longer comparable to them.
        indications = pd.read_csv(indications_path, sep="\t")
        with_insight = indications.loc[
            indications["genetic_insight"] != "none",
            "indication_mesh_id",
        ]
        pp = pp[pp["indication_mesh_id"].isin(set(with_insight))]
        log.info(f"After genetic-insight indication filter: {len(pp):,} programmes")

    pp = pp.assign(phase=pp["ccatnum"].astype(int) - 1)

    # Nothing concluded at the phase a still-running programme currently sits in, so its
    # highest concluded phase is the one below. Reaching a phase concludes every transition
    # beneath it, so running programmes still score successes there. Launched is terminal:
    # the outcome is known, and every launched programme carries it as its active category.
    still_running = pp["acat"].notna() & (pp["acat"] != "Launched")
    pp = pp.assign(phase_concluded=pp["phase"] - still_running)

    pairs = pp.merge(
        load_mesh_to_efo(mesh_to_efo_path),
        left_on="indication_mesh_id",
        right_on="mesh_id",
        how="inner",
    )
    log.info(f"After MeSH to EFO mapping: {len(pairs):,} programme-disease rows")

    # A pair is undetermined only when no programme concluded at its highest phase, so a
    # concluded failure there outranks a sibling programme still running.
    result = pairs.groupby(["gene", "mapped_curie"], as_index=False).agg(
        max_phase=("phase", "max"),
        max_phase_concluded=("phase_concluded", "max"),
    )
    result = result.rename(columns={"mapped_curie": "efo_id"})
    result["is_ongoing"] = result["max_phase_concluded"] < result["max_phase"]
    result = result[["gene", "efo_id", "max_phase", "is_ongoing"]]

    ongoing = int(result["is_ongoing"].sum())
    share = f" ({100 * ongoing / len(result):.1f}%)" if len(result) else ""
    log.info(
        f"Result: {len(result):,} gene-indication pairs, {result['gene'].nunique():,} genes, "
        f"undetermined at highest phase: {ongoing:,}{share}",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def parse_associations(
    associations_path: Path,
    mesh_to_efo_path: Path,
    output_path: Path,
) -> pd.DataFrame:
    """Parse Minikel genetic associations to gene-disease pairs."""
    log.info(f"Loading associations from {associations_path}")
    assoc = pd.read_csv(associations_path, sep="\t", low_memory=False)
    log.info(f"Loaded {len(assoc):,} associations")

    assoc = assoc[assoc["gene"].notna() & assoc["mesh_id"].notna()]

    assoc = assoc[~assoc["source"].isin(SOMATIC_SOURCES)]
    log.info(f"After dropping {sorted(SOMATIC_SOURCES)}: {len(assoc):,}")

    assoc = assoc[assoc["l2g_share"].isna() | (assoc["l2g_share"] >= MIN_L2G_SHARE)]
    log.info(f"After l2g_share >= {MIN_L2G_SHARE}: {len(assoc):,}")

    pairs = assoc.merge(
        load_mesh_to_efo(mesh_to_efo_path),
        on="mesh_id",
        how="inner",
    )
    log.info(f"After MeSH to EFO mapping: {len(pairs):,} association-disease rows")

    result = pairs.groupby(["gene", "mapped_curie"], as_index=False).agg(
        source=("source", lambda x: "&".join(sorted(set(x)))),
    )
    result = result.rename(columns={"mapped_curie": "efo_id"})

    log.info(
        f"Result: {len(result):,} gene-disease pairs, {result['gene'].nunique():,} genes",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_parquet(output_path, index=False)
    log.info(f"Saved to {output_path}")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--pharmaprojects", type=Path, help="Pharmaprojects TSV (pp.tsv)")
    parser.add_argument("--associations", type=Path, help="Associations TSV (assoc.tsv.gz)")
    parser.add_argument(
        "--indications",
        type=Path,
        help="Indication TSV (indic.tsv); restricts to indications with genetic insight",
    )
    parser.add_argument("--mesh-to-efo", type=Path, help="MeSH to EFO mapping TSV")
    parser.add_argument("--output-dir", type=Path, help="Output directory")
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    config_path = args.config or Path(__file__).parent.parent.parent / "configs" / "parsing.yaml"
    if config_path.exists():
        with config_path.open() as f:
            cfg = yaml.safe_load(f)
        ct_inputs = cfg.get("clinical_trials", {})
        ge_inputs = cfg.get("genetic_evidence", {})
        mappings = cfg.get("mappings", {})
        output_dir = Path(cfg.get("output_dir", "."))
    else:
        ct_inputs, ge_inputs, mappings, output_dir = {}, {}, {}, Path()

    mesh_to_efo_path = args.mesh_to_efo or Path(mappings.get("mesh_to_efo", ""))
    out_dir = args.output_dir or output_dir / "minikel"

    parse_clinical_trials(
        pharmaprojects_path=args.pharmaprojects
        or Path(ct_inputs.get("minikel_pharmaprojects", "")),
        mesh_to_efo_path=mesh_to_efo_path,
        output_path=out_dir / "minikel_clinical_trials.parquet",
        indications_path=args.indications or _optional(ct_inputs.get("minikel_indications")),
    )
    parse_associations(
        associations_path=args.associations or Path(ge_inputs.get("minikel_associations", "")),
        mesh_to_efo_path=mesh_to_efo_path,
        output_path=out_dir / "minikel_genetic_evidence.parquet",
    )


if __name__ == "__main__":
    main()
