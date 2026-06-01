"""Configuration schema with defaults."""

from dataclasses import dataclass, field
from pathlib import Path

# Default phase transitions (single source of truth)
DEFAULT_PHASE_TRANSITIONS: tuple[tuple[int, int], ...] = ((1, 2), (2, 3), (3, 4), (1, 4))


@dataclass
class DataPaths:
    """Paths to input data files."""

    clinical_trials: Path | None = None
    genetic_evidence: Path | None = None  # default targets
    baseline_evidence: Path | None = None  # for prioritized mode
    efo_similarity_lookup: Path | None = None
    gene_universe: Path | None = None


@dataclass
class Output:
    """Output configuration."""

    dir: Path | None = None
    save_trials: bool = False
    save_matched_pairs: bool = False


@dataclass
class Thresholds:
    """Validation thresholds with defaults."""

    semantic_similarity: float = 0.8


@dataclass
class Config:
    """Main configuration."""

    data: DataPaths = field(default_factory=DataPaths)
    output: Output = field(default_factory=Output)
    thresholds: Thresholds = field(default_factory=Thresholds)
    phase_transitions: list[list[int]] = field(
        default_factory=lambda: [list(t) for t in DEFAULT_PHASE_TRANSITIONS],
    )


# Phase names - industry standard (0=Preclinical supported)
PHASE_NAMES = {0: "Preclinical", 1: "Phase I", 2: "Phase II", 3: "Phase III", 4: "Approved"}


def phase_label(from_phase: int, to_phase: int) -> str:
    """Generate label for phase transition."""
    from_name = PHASE_NAMES.get(from_phase, str(from_phase))
    to_name = PHASE_NAMES.get(to_phase, str(to_phase))
    # Short forms: Preclinical→Pre, Phase I→I, Approved stays
    from_short = from_name.replace("Phase ", "").replace("Preclinical", "Pre")
    to_short = to_name.replace("Phase ", "").replace("Preclinical", "Pre")
    return f"{from_short}→{to_short}"
