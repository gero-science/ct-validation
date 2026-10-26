"""EFO ontology terms, and remapping of obsolete ids through `replaced_by`.

Every source is coded against whichever EFO release it was built with, while the similarity
lookup is built from one. An id that is not a live term of that release can never be matched
(matching joins through the lookup, diagonal included).

Obsolete terms carry no `is_a`, so they cannot be placed in the lookup; the link from an old
term to its successor is `replaced_by`. `load_efo` follows it transitively, and `remap_efo_ids`
rewrites a source's ids before the source groups by them. Of several `replaced_by`, only the
live ones count. Obsolete terms with no single live successor (only `consider`, a split into
several live terms, nothing, a cycle, or a chain that ends at a non-live term) are dropped: they
cannot be matched.
"""

import logging
import re
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

log = logging.getLogger(__name__)

# efo:EFO_0000094 and IRIs ending in PREFIX_local (EFO, OBO PURLs, ORDO); CURIEs pass through
_PREFIXED = re.compile(r"^(?:.*/|efo:)([A-Z][A-Za-z]*)_([A-Za-z0-9]+)$")


def to_curie(term: str) -> str:
    """Normalise an EFO .obo id or IRI to the colon CURIE the sources use (EFO:0000094)."""
    m = _PREFIXED.match(term)
    return f"{m.group(1)}:{m.group(2)}" if m else term


@dataclass(frozen=True)
class EfoOntology:
    """The terms of one EFO release, keyed by colon CURIE."""

    version: str | None
    live: frozenset[str]
    labels: Mapping[str, str]
    # Obsolete id -> its live successor
    replacements: Mapping[str, str]
    # Obsolete ids with no live successor
    unresolved: frozenset[str]

    def remap(self, ids: pd.Series) -> pd.Series:
        """Replace obsolete ids by their successor and unresolvable ones by NA.

        Live ids, ids this release does not know and NA pass through unchanged.
        """
        replaced = ids.isin(self.replacements.keys())
        dropped = ids.isin(self.unresolved)
        return ids.where(~replaced, ids.map(self.replacements)).mask(dropped)


def _parse_obo(path: Path) -> tuple[str | None, list[dict[str, list[str]]]]:
    """Return the data-version and the [Term] stanzas as {tag: [values]}."""
    version, terms, current = None, [], None
    with path.open() as f:
        for raw in f:
            line = raw.rstrip("\n")
            if line.startswith("["):
                current = {} if line == "[Term]" else None
                if current is not None:
                    terms.append(current)
            elif current is not None and ": " in line:
                tag, value = line.split(": ", 1)
                current.setdefault(tag, []).append(value)
            elif version is None and line.startswith("data-version: "):
                version = line.removeprefix("data-version: ")
    if version and (m := re.search(r"/v(\d+\.\d+\.\d+)/", version)):
        version = m.group(1)
    return version, terms


def load_efo(path: str | Path) -> EfoOntology:
    """Load an EFO .obo: its live terms and a resolved obsolete -> successor map."""
    version, terms = _parse_obo(Path(path))

    live, labels, candidates = set(), {}, {}
    for term in terms:
        curie = to_curie(term["id"][0])
        if term.get("is_obsolete") == ["true"]:
            # The IRI and CURIE spellings of one successor count once
            candidates[curie] = {to_curie(t) for t in term.get("replaced_by", [])}
        else:
            live.add(curie)
            if "name" in term:
                labels[curie] = term["name"][0]
    obsolete = set(candidates)

    replaced_by = {}
    for curie, successors in candidates.items():
        # Of several successors only the live ones count; more than one left is a split
        if len(successors) > 1:
            successors = successors & live
        if len(successors) == 1:
            replaced_by[curie] = next(iter(successors))

    replacements = {}
    for curie in obsolete:
        # Follow the chain until it leaves the obsolete terms, giving up on a cycle
        seen, target = {curie}, replaced_by.get(curie)
        while target in replaced_by and target not in seen:
            seen.add(target)
            target = replaced_by[target]
        if target in live:
            replacements[curie] = target

    log.info(
        f"EFO {version}: {len(live):,} live terms, {len(obsolete):,} obsolete, "
        f"{len(replacements):,} with a live successor",
    )
    return EfoOntology(
        version=version,
        live=frozenset(live),
        labels=labels,
        replacements=replacements,
        unresolved=frozenset(obsolete - replacements.keys()),
    )


def remap_efo_ids(
    df: pd.DataFrame,
    efo: EfoOntology,
    column: str = "efo_id",
    *,
    label_column: str | None = None,
) -> pd.DataFrame:
    """Rewrite obsolete ids in `column` to their successor and drop unresolvable rows.

    Run it before the source groups by disease: two old ids replaced by one term collapse
    into the same pair there. `label_column`, if given, takes the successor's name.
    """
    ids = df[column]
    replaced = ids.isin(efo.replacements.keys())
    dropped = ids.isin(efo.unresolved)
    unknown = ids.notna() & ~ids.isin(efo.live) & ~replaced & ~dropped

    out = df.copy()
    out[column] = efo.remap(ids)
    if label_column is not None:
        out.loc[replaced, label_column] = out.loc[replaced, column].map(efo.labels)
    out = out[~dropped]

    log.info(
        f"EFO {efo.version} remap of {column}: {int(replaced.sum()):,} rows "
        f"({ids[replaced].nunique():,} ids) replaced, {int(dropped.sum()):,} rows "
        f"({ids[dropped].nunique():,} ids) dropped as obsolete with no live successor, "
        f"{int(unknown.sum()):,} rows ({ids[unknown].nunique():,} ids) not in this release",
    )
    return out


def efo_for_remap(mappings: Mapping) -> EfoOntology | None:
    """The EFO to remap against from parsing.yaml `mappings`, or None when remapping is off."""
    if not mappings.get("efo_remap_obsolete"):
        return None
    return load_efo(mappings["efo_ontology"])
