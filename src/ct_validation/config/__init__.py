"""Configuration module."""

from ct_validation.config.loader import load_config
from ct_validation.config.schema import (
    PHASE_NAMES,
    Config,
    DataPaths,
    Output,
    Thresholds,
    phase_label,
)

__all__ = [
    "PHASE_NAMES",
    "Config",
    "DataPaths",
    "Output",
    "Thresholds",
    "load_config",
    "phase_label",
]
