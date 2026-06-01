"""MCP server exposing ct-validation tools for AI agent integration.

Usage:
    python -m ct_validation.mcp_server
    # or via the entry point:
    ct-validation-mcp
"""

import json

import pandas as pd
from mcp.server.fastmcp import FastMCP

from ct_validation import get_expanded_disease_set, validate

mcp = FastMCP(
    "ct-validation",
    instructions=(
        "Clinical trial enrichment validation tools. "
        "Calculates risk ratios for gene-indication pairs progressing through "
        "clinical trial phases, with or without semantic similarity matching. "
        "Accepts file paths (parquet) for large datasets or inline lists for small queries."
    ),
)


def _to_input(
    data: str | list[dict] | None,
    label: str,
) -> pd.DataFrame | str | None:
    """Convert MCP input to DataFrame or file path.

    Accepts:
    - A file path string (ending in .parquet, .csv, .tsv, or containing '/')
    - A JSON string representing a list of dicts
    - A list of dicts directly
    - None
    """
    if data is None:
        return None
    if isinstance(data, list):
        if not data:
            raise ValueError(f"{label}: empty list")
        return pd.DataFrame(data)
    if isinstance(data, str):
        stripped = data.strip()
        if stripped.startswith("["):
            rows = json.loads(stripped)
            if not rows:
                raise ValueError(f"{label}: empty list")
            return pd.DataFrame(rows)
        return stripped
    raise TypeError(f"{label}: expected file path string, JSON array, or list of dicts")


@mcp.tool()
def ct_validate(
    clinical_trials: str | list[dict],
    targets: str | list[dict],
    similarity_lookup: str | list[dict] | None = None,
    baseline_evidence: str | list[dict] | None = None,
    similarity_threshold: float | None = None,
    phase_transitions: list[list[int]] | None = None,
) -> list[dict]:
    """Calculate clinical trial phase transition enrichment for gene-indication pairs.

    Computes risk ratios (RR) and odds ratios (OR) with 95% CI for each phase
    transition, comparing progression rates of gene-indication pairs WITH vs
    WITHOUT genetic evidence.

    Args:
        clinical_trials: Gene-indication pairs with max clinical trial phase.
            File path to parquet, or list of {"gene": str, "efo_id": str, "max_phase": int}.
        targets: Genetic evidence for gene-indication pairs.
            File path to parquet, or list of {"gene": str, "efo_id": str}.
        similarity_lookup: EFO semantic similarity pairs (optional).
            File path to parquet, or list of {"efo_id_1": str, "efo_id_2": str, "similarity": float}.
            If None, performs exact matching on (gene, efo_id) without expansion.
        baseline_evidence: Baseline genetic evidence for prioritized mode (optional).
            Same format as targets. When provided, the "no evidence" group excludes
            pairs with baseline support.
        similarity_threshold: Minimum semantic similarity (default 0.8).
            Only used when similarity_lookup is provided.
        phase_transitions: Phase transitions to analyze as [[from, to], ...].
            Default: [[1,2], [2,3], [3,4], [1,4]].

    Returns:
        List of dicts with enrichment results per phase transition:
        - phase_from, phase_to, phase_label
        - n_yes, n_no, x_yes, x_no, rate_yes, rate_no
        - rr, rr_ci_lower, rr_ci_upper (risk ratio with 95% CI)
        - or, or_ci_lower, or_ci_upper (odds ratio with 95% CI)
        - p_value (two-sided Fisher's exact test)
    """
    ct_input = _to_input(clinical_trials, "clinical_trials")
    targets_input = _to_input(targets, "targets")
    sim_input = _to_input(similarity_lookup, "similarity_lookup")
    baseline_input = _to_input(baseline_evidence, "baseline_evidence")

    parsed_transitions = None
    if phase_transitions is not None:
        parsed_transitions = [tuple(t) for t in phase_transitions]

    result = validate(
        clinical_trials=ct_input,
        targets=targets_input,
        similarity_lookup=sim_input,
        baseline_evidence=baseline_input,
        similarity_threshold=similarity_threshold,
        phase_transitions=parsed_transitions,
    )

    return result.to_dict(orient="records")


@mcp.tool()
def expand_disease_set(
    efo_ids: list[str],
    similarity_lookup: str | list[dict],
    similarity_threshold: float = 0.8,
) -> list[str]:
    """Expand a set of EFO disease IDs using semantic similarity.

    Returns the input diseases plus similar diseases from the similarity matrix
    that meet the threshold. Useful for filtering clinical trials to a disease
    neighborhood before running validation.

    Args:
        efo_ids: List of EFO disease identifiers (e.g., ["EFO:0000270", "EFO:0000384"]).
        similarity_lookup: EFO semantic similarity pairs.
            File path to parquet, or list of {"efo_id_1": str, "efo_id_2": str, "similarity": float}.
        similarity_threshold: Minimum similarity to include (default 0.8).

    Returns:
        List of EFO IDs: input diseases + semantically similar diseases.
    """
    sim_input = _to_input(similarity_lookup, "similarity_lookup")

    expanded = get_expanded_disease_set(
        efo_ids=set(efo_ids),
        similarity_pairs=sim_input,
        similarity_threshold=similarity_threshold,
    )

    return sorted(expanded)


def main() -> None:
    """Run the MCP server."""
    mcp.run()


if __name__ == "__main__":
    main()
