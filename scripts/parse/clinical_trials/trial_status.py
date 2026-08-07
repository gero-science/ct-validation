#!/usr/bin/env python3
"""Resolve clinical-trial record status into reached/concluded facts.

A record carries two independent facts: the programme reached the phase, and
whether its outcome there is determined. Only Withdrawn fails the first — those
trials enrolled nobody, so the record cannot show the phase was reached.
Dropping any other status would demote pairs into fabricated failures one phase
down.

Source vocabularies differ but the classification must not, so both live here.
"""

import logging

import pandas as pd

log = logging.getLogger(__name__)

# Status -> outcome at this phase is determined.
OPENTARGETS_CONCLUDED = {
    "Completed": True,
    "Terminated": True,
    # ClinicalTrials.gov derives this only after the completion date has passed,
    # so these went dark rather than stalled
    "Unknown status": True,
    None: True,  # approval records (DailyMed/ATC/EMA/FDA), not trials
    "Recruiting": False,
    "Active, not recruiting": False,
    "Not yet recruiting": False,
    "Enrolling by invitation": False,
    "Suspended": False,  # halted, but participants enrolled; outcome still open
}
OPENTARGETS_DROPPED = {"Withdrawn"}

TRIALPANORAMA_CONCLUDED = {
    "COMPLETED": True,
    "TERMINATED": True,
    "UNKNOWN": True,
    "RECRUITING": False,
    "ACTIVE_NOT_RECRUITING": False,
    "NOT_YET_RECRUITING": False,
    "ENROLLING_BY_INVITATION": False,
    "SUSPENDED": False,
}
TRIALPANORAMA_DROPPED = {"WITHDRAWN"}

# ClinicalTrials.gov's Unknown label, in both source vocabularies. Read as concluded
# above; aggregate.py's --unknown-ongoing is the sensitivity arm for that choice.
UNKNOWN_STATUSES = {"Unknown status", "UNKNOWN"}


def resolve(
    status: pd.Series,
    concluded: dict,
    dropped: set,
) -> tuple[pd.Series, pd.Series]:
    """Classify trial records by status.

    Args:
        concluded: Status -> whether the outcome at this phase is determined.
        dropped: Statuses that do not establish the phase was reached.

    Returns:
        (keep, is_concluded), aligned to `status`. `is_concluded` is only
        meaningful where `keep`.

    Raises:
        ValueError: on a status in neither mapping — a new status must be
            classified deliberately rather than silently censoring rows.
    """
    present = set(status.dropna().unique())
    if status.isna().any():
        present.add(None)
    if unmapped := present - set(concluded) - dropped:
        raise ValueError(f"unclassified trial status: {sorted(map(str, unmapped))}")

    keep = ~status.isin(dropped)
    is_concluded = status.map(concluded).astype("boolean")
    if None in concluded:
        # .map() matches a literal None key but not NaN, so nulls need masking
        is_concluded = is_concluded.mask(status.isna(), concluded[None])
    return keep, is_concluded


def apply(df: pd.DataFrame, concluded: dict, dropped: set) -> pd.DataFrame:
    """Drop records that do not establish a phase, and add `is_concluded`."""
    keep, is_concluded = resolve(df["status"], concluded, dropped)
    log.info(f"Dropping {(~keep).sum():,} records with status in {sorted(dropped)}")
    records = df[keep].copy()
    records["is_concluded"] = is_concluded[keep].astype(bool)
    return records
