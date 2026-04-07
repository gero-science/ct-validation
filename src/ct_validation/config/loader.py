"""Config loader from YAML files."""

from pathlib import Path

import yaml

from ct_validation.config.schema import Config, DataPaths, Output, Thresholds


def load_config(path: str | Path) -> Config:
    """Load config from YAML file. Missing keys use dataclass defaults."""
    with Path(path).open() as f:
        raw = yaml.safe_load(f) or {}

    # Parse data paths - convert strings to Path objects
    data_raw = raw.get("data", {})
    data_kwargs = {k: Path(v) for k, v in data_raw.items() if v}
    data = DataPaths(**data_kwargs)

    # Parse output config
    output_raw = raw.get("output", {})
    output_kwargs = {}
    if output_raw.get("dir"):
        output_kwargs["dir"] = Path(output_raw["dir"])
    if "save_trials" in output_raw:
        output_kwargs["save_trials"] = output_raw["save_trials"]
    output = Output(**output_kwargs)

    # Parse thresholds - dataclass handles defaults
    thresh_raw = raw.get("thresholds", {})
    thresholds = Thresholds(**thresh_raw)

    # Parse phase transitions - dataclass handles default
    phase_transitions = raw.get("phase_transitions")

    return Config(
        data=data,
        output=output,
        thresholds=thresholds,
        **({"phase_transitions": phase_transitions} if phase_transitions else {}),
    )
