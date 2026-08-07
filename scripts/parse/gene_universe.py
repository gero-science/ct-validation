#!/usr/bin/env python3
"""Write the protein-coding gene universe from the OpenTargets target table.

One approved symbol per line, matching the gene symbols the parsers emit.
"""

import argparse
import logging
from pathlib import Path

import pandas as pd
import yaml

log = logging.getLogger(__name__)


def build_gene_universe(target_path: Path, output_path: Path) -> list[str]:
    """Write sorted protein-coding gene symbols from the OpenTargets target table."""
    log.info(f"Loading targets from {target_path}")
    targets = pd.read_parquet(target_path, columns=["approvedSymbol", "biotype"])
    log.info(f"Loaded {len(targets):,} targets")

    protein_coding = targets[targets["biotype"] == "protein_coding"]
    symbols = sorted(set(protein_coding["approvedSymbol"].dropna()))
    # approvedSymbol falls back to the Ensembl ID where a gene has no HGNC symbol; the
    # genetic-evidence parser names those genes the same way, so they are kept.
    n_unnamed = sum(s.startswith("ENSG") for s in symbols)
    log.info(
        f"Kept {len(protein_coding):,} protein-coding targets, {len(symbols):,} symbols "
        f"({n_unnamed:,} unnamed, carrying an Ensembl ID)",
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text("\n".join(symbols) + "\n")
    log.info(f"Saved to {output_path}")
    return symbols


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, help="Config YAML")
    parser.add_argument("--target", type=Path, help="OpenTargets target parquet directory")
    parser.add_argument("--output", type=Path, help="Output text path")
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
        ge_inputs = cfg.get("genetic_evidence", {})
        mappings = cfg.get("mappings", {})
    else:
        ge_inputs, mappings = {}, {}

    build_gene_universe(
        target_path=args.target or Path(ge_inputs.get("opentargets_target", "")),
        output_path=args.output or Path(mappings.get("gene_universe", "")),
    )


if __name__ == "__main__":
    main()
