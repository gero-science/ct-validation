#!/usr/bin/env python3
"""Run parsing scripts based on config.

Usage:
    python scripts/parsing/run_parsing.py                    # Run all enabled
    python scripts/parsing/run_parsing.py --only genetic_evidence
    python scripts/parsing/run_parsing.py --only clinical_trials --skip-aggregate
    python scripts/parsing/run_parsing.py --only clinical_trials --sources stitch

"""

import argparse
import logging
import subprocess
import sys
from pathlib import Path

import yaml

log = logging.getLogger(__name__)

SCRIPTS_DIR = Path(__file__).parent

GENETIC_EVIDENCE_SCRIPTS = {
    "gwas_catalog": SCRIPTS_DIR / "genetic_evidence" / "gwas_catalog.py",
    "clinvar": SCRIPTS_DIR / "genetic_evidence" / "clinvar.py",
    "opentargets": SCRIPTS_DIR / "genetic_evidence" / "opentargets.py",
    "omim": SCRIPTS_DIR / "genetic_evidence" / "omim.py",
    "genebass": SCRIPTS_DIR / "genetic_evidence" / "genebass.py",
}

CLINICAL_TRIALS_SCRIPTS = {
    "stitch": SCRIPTS_DIR / "clinical_trials" / "stitch.py",
    "dgidb": SCRIPTS_DIR / "clinical_trials" / "dgidb.py",
    "trialpanorama": SCRIPTS_DIR / "clinical_trials" / "trialpanorama.py",
    "chembl": SCRIPTS_DIR / "clinical_trials" / "chembl.py",
    "opentargets": SCRIPTS_DIR / "clinical_trials" / "opentargets.py",
}

AGGREGATE_SCRIPTS = {
    "genetic_evidence": SCRIPTS_DIR / "genetic_evidence" / "aggregate.py",
    "clinical_trials": SCRIPTS_DIR / "clinical_trials" / "aggregate.py",
}


def run_script(
    script_path: Path,
    *,
    verbose: bool = False,
    extra_args: list[str] | None = None,
) -> bool:
    """Run a script and return success status."""
    cmd = [sys.executable, str(script_path)]
    if verbose:
        cmd.append("-v")
    if extra_args:
        cmd.extend(extra_args)

    log.info(
        f"Running: {script_path.name}"
        + (f" {' '.join(extra_args)}" if extra_args else ""),
    )
    result = subprocess.run(cmd, check=False)

    if result.returncode != 0:
        log.error(f"FAILED: {script_path.name} (exit code {result.returncode})")
        return False
    return True


def run_domain(
    domain: str,
    scripts: dict[str, Path],
    sources: list[str] | None,
    *,
    aggregate: bool,
    verbose: bool,
    script_args: dict[str, list[str]] | None = None,
) -> bool:
    """Run scripts for a domain."""
    log.info(f"=== {domain.upper()} ===")
    script_args = script_args or {}

    # Run source parsers
    sources_to_run = sources or list(scripts.keys())
    for source in sources_to_run:
        if source not in scripts:
            log.warning(f"Unknown source: {source}")
            continue
        if not scripts[source].exists():
            log.warning(f"Script not found: {scripts[source]}")
            continue
        extra = script_args.get(source)
        if not run_script(scripts[source], verbose=verbose, extra_args=extra):
            return False

    # Run aggregation
    if aggregate and domain in AGGREGATE_SCRIPTS:
        if not run_script(AGGREGATE_SCRIPTS[domain], verbose=verbose):
            return False
    return True


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument(
        "--only",
        choices=["genetic_evidence", "clinical_trials"],
        help="Run only this domain",
    )
    parser.add_argument("--sources", type=str, help="Comma-separated sources to run")
    parser.add_argument(
        "--skip-aggregate", action="store_true", help="Skip aggregation step",
    )
    parser.add_argument("-v", "--verbose", action="store_true")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    # Load config
    config_path = args.config or SCRIPTS_DIR.parent.parent / "configs" / "parsing.yaml"
    if config_path.exists():
        with config_path.open() as f:
            cfg = yaml.safe_load(f)
        run_cfg = cfg.get("run", {})
    else:
        log.warning(f"Config not found: {config_path}")
        run_cfg = {}

    # Parse CLI sources
    cli_sources = args.sources.split(",") if args.sources else None

    # Run genetic evidence
    ge_cfg = run_cfg.get("genetic_evidence", {})
    if (args.only is None or args.only == "genetic_evidence") and ge_cfg.get(
        "enabled", True,
    ):
        sources = cli_sources or ge_cfg.get("sources")
        if not run_domain(
            "genetic_evidence",
            GENETIC_EVIDENCE_SCRIPTS,
            sources,
            aggregate=not args.skip_aggregate,
            verbose=args.verbose,
        ):
            sys.exit(1)

    # Run clinical trials
    ct_cfg = run_cfg.get("clinical_trials", {})
    if (args.only is None or args.only == "clinical_trials") and ct_cfg.get(
        "enabled", True,
    ):
        sources = cli_sources or ct_cfg.get("sources")
        if not run_domain(
            "clinical_trials",
            CLINICAL_TRIALS_SCRIPTS,
            sources,
            aggregate=not args.skip_aggregate,
            verbose=args.verbose,
        ):
            sys.exit(1)

    log.info("Done.")


if __name__ == "__main__":
    main()
