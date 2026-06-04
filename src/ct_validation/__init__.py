"""Clinical trial validation package.

Main entry point:
- validate(): Run validation pipeline (with config and/or explicit args)
"""

from ct_validation.api import validate
from ct_validation.config import load_config
from ct_validation.config.schema import PHASE_NAMES, Config, phase_label
from ct_validation.plotting import forest_plot
from ct_validation.validation import (
    EnrichmentResult,
    calculate_enrichment,
    create_matched_pairs_df,
    create_matched_pairs_set,
    get_expanded_disease_set,
    katz_ci_risk_ratio,
    woolf_ci_odds_ratio,
)

__version__ = "0.2.0"

__all__ = [
    "PHASE_NAMES",
    "Config",
    "EnrichmentResult",
    "calculate_enrichment",
    "create_matched_pairs_df",
    "create_matched_pairs_set",
    "forest_plot",
    "get_expanded_disease_set",
    "katz_ci_risk_ratio",
    "load_config",
    "phase_label",
    "validate",  # Main API
    "woolf_ci_odds_ratio",
]
