"""Command-line interface for ct-validation."""

from enum import Enum
from pathlib import Path
from typing import Annotated

import typer

from ct_validation import validate
from ct_validation.config import load_config


class OutputFormat(str, Enum):
    parquet = "parquet"
    tsv = "tsv"
    csv = "csv"


cli = typer.Typer(add_completion=False)


@cli.command()
def main(
    config_path: Annotated[
        Path | None,
        typer.Option("--config", "-c", exists=True, help="Path to YAML config file."),
    ] = None,
    clinical_trials: Annotated[
        Path | None,
        typer.Option(exists=True, help="Path to clinical trials parquet."),
    ] = None,
    targets: Annotated[
        list[Path] | None,
        typer.Option(
            exists=True,
            help="Path to genetic evidence (targets) parquet. Repeat for batch.",
        ),
    ] = None,
    similarity_lookup: Annotated[
        Path | None,
        typer.Option(exists=True, help="Path to EFO similarity lookup parquet."),
    ] = None,
    baseline_evidence: Annotated[
        Path | None,
        typer.Option(exists=True, help="Path to baseline evidence parquet (for prioritized mode)."),
    ] = None,
    gene_universe: Annotated[
        Path | None,
        typer.Option(exists=True, help="Path to gene universe file (one gene per line)."),
    ] = None,
    similarity_threshold: Annotated[
        float | None,
        typer.Option(help="Minimum semantic similarity threshold (default: 0.8)."),
    ] = None,
    phase_transitions: Annotated[
        str | None,
        typer.Option(
            help="Phase transitions to analyze, e.g. '1,2;2,3;3,4;1,4' (default: 1,2;2,3;3,4;1,4).",
        ),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option("-o", "--output-dir", help="Output directory for results."),
    ] = None,
    save_trials: Annotated[  # noqa: FBT002
        bool,
        typer.Option(
            "--save-trials",
            help="Save annotated trials DataFrame alongside enrichment results.",
        ),
    ] = False,
    output_format: Annotated[
        OutputFormat,
        typer.Option("--format", help="Output format."),
    ] = OutputFormat.parquet,
) -> None:
    """Run clinical trial enrichment validation.

    Can be called with config only, args only, or both (args override config).

    Input Schemas:
        clinical_trials:   gene, efo_id, max_phase
        targets:           gene, efo_id
        similarity_lookup: efo_id_1, efo_id_2, similarity
        baseline_evidence: gene, efo_id (optional)

    Output Schema:
        phase_from, phase_to    Phase transition
        n_yes, n_no             Pairs at phase_from+ (with/without evidence)
        x_yes, x_no             Pairs reaching phase_to+
        rate_yes, rate_no       Progression rates
        rr, rr_ci_lower/upper   Risk ratio with 95% CI
        or, or_ci_lower/upper   Odds ratio with 95% CI

    Examples:
        ct-validation --config config.yaml
        ct-validation --clinical-trials ct.parquet --targets ge.parquet --similarity-lookup sim.parquet -o results/
        ct-validation --config config.yaml --similarity-threshold 0.9 -o custom_output/
    """
    # Load config if provided
    config = None
    if config_path:
        config = load_config(config_path)

    # Resolve output directory
    if output_dir is None:
        output_dir = config.output.dir if config and config.output.dir else Path()

    # Resolve save_trials from config if not set via flag
    if not save_trials and config:
        save_trials = config.output.save_trials

    output_dir.mkdir(parents=True, exist_ok=True)

    # Parse phase transitions
    parsed_transitions = None
    if phase_transitions:
        parsed_transitions = [
            tuple(int(x) for x in pair.split(",")) for pair in phase_transitions.split(";")
        ]

    # Normalize targets: single Path -> pass as-is, multiple -> pass as list
    targets_arg: Path | list[Path] | None = None
    if targets is not None:
        targets_arg = targets[0] if len(targets) == 1 else targets
    batch = isinstance(targets_arg, list)

    # Run validation
    result = validate(
        config=config,
        clinical_trials=clinical_trials,
        targets=targets_arg,
        similarity_lookup=similarity_lookup,
        baseline_evidence=baseline_evidence,
        gene_universe=gene_universe,
        similarity_threshold=similarity_threshold,
        phase_transitions=parsed_transitions,
        return_trials=save_trials,
    )

    # Normalize to list for uniform handling
    results = result if batch else [result]
    labels = [t.stem for t in targets] if batch else [None]

    for r, label in zip(results, labels):
        suffix = f"_{label}" if label else ""

        if save_trials:
            enrichment_df, trials_df = r
        else:
            enrichment_df = r

        # Save results
        enrichment_path = output_dir / f"enrichment_results{suffix}.{output_format.value}"
        _save_df(enrichment_df, enrichment_path, output_format)
        print(f"Saved enrichment results to {enrichment_path}")

        if save_trials:
            trials_path = output_dir / f"annotated_trials{suffix}.{output_format.value}"
            _save_df(trials_df, trials_path, output_format)
            print(f"Saved annotated trials to {trials_path}")

        # Print summary
        header = f"\nEnrichment Results ({label}):" if label else "\nEnrichment Results:"
        print(header)
        print("-" * 60)
        print(f"{'Phase':<12} {'n_yes':>8} {'n_no':>10} {'RR':>8} {'95% CI':<20}")
        print("-" * 60)
        for _, row in enrichment_df.iterrows():
            ci = f"[{row['rr_ci_lower']:.3f}, {row['rr_ci_upper']:.3f}]"
            print(
                f"{row['phase_label']:<12} {row['n_yes']:>8} {row['n_no']:>10} "
                f"{row['rr']:>8.3f} {ci:<20}",
            )


def _save_df(df, path: Path, fmt: OutputFormat) -> None:
    """Save DataFrame in specified format."""
    if fmt == OutputFormat.parquet:
        df.to_parquet(path, index=False)
    elif fmt == OutputFormat.tsv:
        df.to_csv(path, sep="\t", index=False)
    elif fmt == OutputFormat.csv:
        df.to_csv(path, index=False)


if __name__ == "__main__":
    cli()
