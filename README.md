# <img width="300" alt="ct-validation" src="https://github.com/user-attachments/assets/6dc30335-d5f3-432c-b09c-e46a11f2b575" />

An open framework for benchmarking target–indication evidence against clinical trial outcomes.

[![PyPI](https://img.shields.io/pypi/v/ct-validation)](https://pypi.org/project/ct-validation/)
[![Python](https://img.shields.io/pypi/pyversions/ct-validation)](https://pypi.org/project/ct-validation/)
[![CI](https://github.com/gero-science/ct-validation/actions/workflows/ci.yml/badge.svg)](https://github.com/gero-science/ct-validation/actions/workflows/ci.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21840865.svg)](https://doi.org/10.5281/zenodo.21840865)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

`ct-validation` tests whether target–indication pairs supported by a given evidence source are
enriched for clinical progression. The evidence may be genetic, omics-derived, computational or
curated — any set of `(gene, disease)` pairs.

```python
import ct_validation as ctv

ctv.validate(
    clinical_trials="data/clinical_trials/aggregated/gene_indication_max_phase.parquet",
    targets="data/genetic_evidence/aggregated/genetic_evidence.parquet",
    similarity_lookup="data/mappings/efo_similarity_lookup_0.5.parquet",
)
```

```text
 phase_label  x_yes  n_yes    x_no     n_no    rr  rr_ci_lower  rr_ci_upper    p_value
        I→II  17541  18544 1934679  2118134  1.04         1.03         1.04   6.5e-63
      II→III  12979  16603 1184685  1785165  1.18         1.17         1.19  1.0e-242
III→Approved   8576  12221  623052  1076226  1.21         1.20         1.23  4.7e-171
  I→Approved   8576  16848  623052  1860161  1.52         1.50         1.54        0.0
```

The final row is the phase I to approval composite: supported pairs progressed **1.52 times as
often** as unsupported ones (95% CI 1.50–1.54).
`rr` is the *relative success* (RS) metric of [Nelson et al.
(2015)](https://doi.org/10.1038/ng.3314) and [Minikel et al.
(2024)](https://doi.org/10.1038/s41586-024-07316-0) — a risk ratio of phase advancement.

## Install

```bash
pip install ct-validation          # Python 3.11+
pip install ct-validation[plot]    # + forest plots
pip install ct-validation[mcp]     # + MCP server for agent workflows
```

## Applications

### Benchmarking a single evidence source

As in the example above. Any table carrying `gene` and `efo_id` columns is accepted; the framework is
agnostic to the provenance of the evidence.

### Comparing evidence sources

```python
results = ctv.validate(
    clinical_trials="ct.parquet",
    targets=["gwas_catalog.parquet", "clinvar.parquet", "genebass.parquet"],
    similarity_lookup="efo_similarity_lookup_0.5.parquet",
)
# one DataFrame per source
```

### Incremental value over a baseline

```python
results = ctv.validate(
    clinical_trials="ct.parquet",
    targets="novel_score.parquet",
    baseline_evidence="established_genetics.parquet",   # baseline-only pairs are excluded
)
```

Pairs supported only by the baseline are excluded, so the estimate reflects the incremental value
of the new source rather than its overlap with established evidence.

### Semantic disease matching

Evidence rarely uses the exact EFO term a trial was registered under. When a similarity lookup is
supplied, a clinical pair `(gene, disease)` is supported if the evidence contains `(gene, disease′)`
whose disease scores at or above the matching threshold (default 0.8; exact matches score 1):

```python
expanded = ctv.get_expanded_disease_set(
    efo_ids={"EFO:0000270", "EFO:0000384"},          # asthma, Crohn's disease
    similarity_pairs="efo_similarity_lookup_0.5.parquet",
    similarity_threshold=0.8,
)
```

### Visualization

```python
results = ctv.validate(...)
ctv.forest_plot(results, metric="rr", title="Phase I → Approved")
```

## Other interfaces

<details>
<summary><b>Command line</b></summary>

```bash
# from a config file
ct-validation --config configs/default.yaml

# or explicit arguments
ct-validation \
    --clinical-trials ct.parquet \
    --targets evidence.parquet \
    --similarity-lookup similarity.parquet \
    -o results/

# multiple sources at once
ct-validation --clinical-trials ct.parquet \
    --targets gwas.parquet --targets clinvar.parquet --targets genebass.parquet \
    -o results/

# also write the annotated trials and the matched-pair audit table
ct-validation --config configs/default.yaml --save-trials --save-matched-pairs -o results/
```

See `configs/default.yaml` for validation settings and `configs/parsing.yaml` for data source
paths. Every config value can be overridden by a CLI argument.

</details>

<details>
<summary><b>MCP server</b> (agent workflows)</summary>

```bash
ct-validation-mcp
```

Two tools:

- `ct_validate` — phase-transition enrichment, from Parquet paths or inline records
- `expand_disease_set` — expand EFO IDs through semantic similarity

</details>

## Inputs and outputs

<details>
<summary><b>Input schemas</b></summary>

| Input               | Columns                                                | Description                                                                                 |
| ------------------- | ------------------------------------------------------ | ------------------------------------------------------------------------------------------- |
| `clinical_trials`   | `gene`, `efo_id`, `max_phase`, `is_ongoing` (optional) | Target–indication pairs with highest phase reached, and whether that phase is still running |
| `targets`           | `gene`, `efo_id`                                       | Gene–disease pairs with supporting evidence                                                 |
| `similarity_lookup` | `efo_id_1`, `efo_id_2`, `similarity`                   | Pairwise EFO similarity, **symmetric with diagonal** (optional)                             |
| `baseline_evidence` | `gene`, `efo_id`                                       | Baseline evidence for prioritized mode (optional)                                           |
| `gene_universe`     | one gene per line (text file)                          | Restrict analysis to these genes (optional)                                                 |

All inputs accept Parquet paths or pandas DataFrames (except `gene_universe`, a text file or a
Python set).

`clinical_trials` must carry **one row per `(gene, efo_id)` pair** — rows are counted
individually, so a repeated pair would inflate the denominators. `validate()` raises on duplicates
rather than collapsing them, since only the caller knows what `max_phase` and `is_ongoing` should
be for a merged row.

The similarity lookup **must be symmetric**: matching probes a single orientation, so every pair
`(a, b)` needs `(b, a)` too, plus a self-similarity diagonal. Check one with
`ctv.check_similarity_lookup(sim)`, or pass `check_similarity=True` to `validate()`; it raises
rather than silently dropping half of each disease neighbourhood.

</details>

<details>
<summary><b>Output schema</b></summary>

| Column                             | Description                                                 |
| ---------------------------------- | ----------------------------------------------------------- |
| `phase_from`, `phase_to`           | Phase transition (e.g. 1→2, 1→4)                            |
| `phase_label`                      | Human-readable label (e.g. `I→Approved`)                    |
| `x_yes`, `x_no`                    | Pairs reaching the target phase (with/without evidence)     |
| `n_yes`, `n_no`                    | Pairs entering the transition, less censored                |
| `rate_yes`, `rate_no`              | Progression rates                                           |
| `rr`, `rr_ci_lower`, `rr_ci_upper` | Relative success (risk ratio) with 95% CI (Katz log method) |
| `or`, `or_ci_lower`, `or_ci_upper` | Odds ratio with 95% CI (Woolf logit method)                 |
| `p_value`                          | Two-sided Fisher's exact test p-value                       |

When either group is empty (`n_yes=0` or `n_no=0`), the ratios, their intervals and `p_value` are
`NaN`. A Haldane–Anscombe correction of 0.5 is applied to all four cells when any one of them is
zero.

</details>

<details>
<summary><b>Exporting annotated trials and matched pairs</b></summary>

```python
enrichment, trials = ctv.validate(..., return_trials=True)
enrichment, matched = ctv.validate(..., return_matched_pairs=True)
enrichment, trials, matched = ctv.validate(..., return_trials=True, return_matched_pairs=True)
```

The match-audit table can also be built on its own:

```python
matched = ctv.create_matched_pairs_df(
    genetic_evidence=..., clinical_trials=..., similarity_pairs=...,
    similarity_threshold=0.8,
)
```

| Column       | Description                                              |
| ------------ | -------------------------------------------------------- |
| `gene`       | Gene symbol                                              |
| `ct_efo_id`  | Disease on the clinical trial row                        |
| `ge_efo_id`  | Supporting evidence disease                              |
| `similarity` | Match score (`1.0` for exact matches)                    |
| `match_type` | `exact` when `ct_efo_id == ge_efo_id`, else `similarity` |

</details>

## How the enrichment is computed

For each transition, pairs that reached at least the starting phase are split into supported and
unsupported groups, and

```
RS = (x_yes / n_yes) / (x_no / n_no)
```

Greater than one means supported pairs progress more often.

<details>
<summary><b>Censoring undetermined outcomes</b></summary>

When `clinical_trials` carries an `is_ongoing` column, a pair whose highest phase is still running
is censored from any transition it has not yet completed — dropped from both `n` and `x` where
`is_ongoing` and `max_phase < phase_to`. Such a pair has not failed the transition; its outcome is
not known yet. Censoring is per transition, so a pair that already reached `phase_to` still counts
as a success.

Without the column every pair counts as concluded, nothing is censored, and `validate()` warns
(`MissingOngoingWarning`) — a legitimate choice, and the only one available when the source
carries no trial status. A missing value reads as not-ongoing; non-boolean values raise.

`is_ongoing` is the caller's to define. The parsers in `scripts/` set it where nothing concluded
at the pair's highest phase **and** that phase is below approval — reaching approval settles the
outcome however many post-marketing trials are still recruiting.

</details>

<details>
<summary><b>Clustered bootstrap intervals</b></summary>

Katz and Woolf assume independent rows, but one gene spans many indications — so their intervals
run narrow when the correlated unit is the gene rather than the pair. Pass `cluster_col` to add
percentile intervals that resample clusters with replacement instead of rows:

```python
results = ctv.validate(..., cluster_col="gene")
```

This adds four columns — `rr_boot_ci_lower`, `rr_boot_ci_upper`, `or_boot_ci_lower`,
`or_boot_ci_upper` — and nothing else: resampling only re-estimates spread, so point estimates,
counts and closed-form intervals are unchanged. `p_value` is unaffected too — it stays a row-level
Fisher's exact test on the observed table, so it does not carry the clustering the intervals do.

`bootstrap_replicates` (default 10,000) and `bootstrap_seed` (default 0) control the resampling.
The percentile bounds carry Monte Carlo error falling as `1/√n`; at the default each bound shifts
by a few percent of the interval's width across seeds. Both are ignored without `cluster_col`.

No continuity correction is applied, unlike the closed-form intervals — a percentile bootstrap
needs neither `log(ratio)` nor its standard error. A replicate that empties a cell keeps its
ratio, so **a bound may be returned as `0` or infinite**; relevant when log-scaling these columns.

Two clusters per arm is the floor for a defined interval, not a threshold at which one is
trustworthy — a percentile cluster bootstrap undercovers badly on a few dozen clusters, so read
bounds from a small cluster count as indicative.

Bounds are `NaN` when an arm holds fewer than two clusters, when the observed table has a zero
cell, or when over 25% of replicates leave a ratio undefined (`0/0`) — past 5% a
`BootstrapReplicateLossWarning` is emitted. The closed-form columns stay populated unless an arm is empty
outright.

</details>

## Data

The `scripts/` directory contains reproducible parsers that convert public databases into the two
input tables, harmonized on EFO.

<details>
<summary><b>What each parser covers</b></summary>

**Genetic evidence** (`scripts/parse/genetic_evidence/`):

- GWAS Catalog — genome-wide significant associations (p < 5e-8)
- ClinVar — pathogenic/likely pathogenic variants at high review confidence (somatic,
  susceptibility and response entries excluded)
- Open Targets — germline genetic evidence streams (direct score ≥ 0.5); somatic sources
  (`eva_somatic`, `intogen`) are opt-in via `--include-somatic`
- Genebass — exome-wide associations (p < 1e-7)

**Clinical trials** (`scripts/parse/clinical_trials/`):

- ChEMBL — gene–drug links (pChEMBL ≥ 7.0); its indications are parsed but excluded from the
  aggregate, since they carry no trial status to censor on
- Open Targets — known drug and indication data
- STITCH — high-confidence activation/inhibition links
- DGIdb — drug–gene interactions
- TrialPanorama — interventional studies

**Ontology** (`scripts/r/`):

- EFO semantic similarity matrix (averaged normalized Resnik + Lin information content)

**Benchmark comparator** (`scripts/parse/minikel.py`):

- Minikel et al. (2024) — Citeline programmes with per-transition censoring, and their genetic
  associations, for the published-vs-ours comparison

Configure paths in `configs/parsing.yaml`, then:

```bash
python scripts/parse/gene_universe.py
python scripts/parse/minikel.py
python scripts/parse/run_parsing.py
```

See [DATA_SOURCES.md](DATA_SOURCES.md) for download links, versions, licences and fetching
instructions.

</details>

<details>
<summary><b>Reproducing the paper</b></summary>

The data deposit carries every input `notebooks/benchmark.py` reads, so nothing has to be fetched:

```bash
tar xzf ct-validation-data.tar.gz          # at the repository root
pip install -r notebooks/requirements.txt  # pinned env: Jupyter + plotting theme + editable install
cd notebooks && python benchmark.py
```

The full run takes about 10 minutes and peaks at 8.5 GB resident memory, so 16 GB of RAM is
recommended.

| Resource                      | Location                                                           |
| ----------------------------- | ------------------------------------------------------------------ |
| Source code                   | https://github.com/gero-science/ct-validation                      |
| Archived code (all versions)  | [10.5281/zenodo.21840865](https://doi.org/10.5281/zenodo.21840865) |
| Benchmark data (all versions) | [10.5281/zenodo.21839216](https://doi.org/10.5281/zenodo.21839216) |
| Package                       | [PyPI](https://pypi.org/project/ct-validation/)                    |

</details>

<details>
<summary><b>Development</b></summary>

CI runs on push/PR to `main` (Python 3.11–3.13, `ruff`, `pytest`):

```bash
uv sync --extra mcp --extra parse
uv run ruff check src tests scripts
uv run pytest tests/ -v
```

</details>

## Citation

> Kostiuk K, Igumnov D, Fedichev P, Feizi A. _ct-validation: an open framework for benchmarking
> target–indication evidence against clinical trial outcomes._ (2026)

Archived release [10.5281/zenodo.21922222](https://doi.org/10.5281/zenodo.21922222) · benchmark
data [10.5281/zenodo.21839217](https://doi.org/10.5281/zenodo.21839217) · machine-readable
metadata in [CITATION.cff](CITATION.cff).

## License

MIT
