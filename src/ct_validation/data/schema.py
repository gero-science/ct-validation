"""Data schemas for validation pipeline."""

from dataclasses import dataclass
from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


@dataclass(frozen=True)
class Schema:
    """Schema definition with required columns."""

    required: tuple[str, ...]

    def validate(self, df: pd.DataFrame, name: str = "data") -> None:
        """Raise ValueError if required columns are missing."""
        missing = [c for c in self.required if c not in df.columns]
        if missing:
            raise ValueError(f"{name} missing required columns: {missing}")

    def validate_parquet(self, path: str | Path, name: str = "data") -> None:
        """Validate parquet schema from metadata only (no data loading)."""
        pf = pq.ParquetFile(path)
        columns = pf.schema_arrow.names
        missing = [c for c in self.required if c not in columns]
        if missing:
            raise ValueError(f"{name} missing required columns: {missing}")


GENETIC_EVIDENCE = Schema(required=("gene", "efo_id"))
CLINICAL_TRIALS = Schema(required=("gene", "efo_id", "max_phase"))
SIMILARITY_LOOKUP = Schema(required=("efo_id_1", "efo_id_2", "similarity"))
