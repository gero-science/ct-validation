# <img width="300" alt="ct-validation" src="https://github.com/user-attachments/assets/6dc30335-d5f3-432c-b09c-e46a11f2b575" />

An open framework for benchmarking gene-indication evidence against clinical trial outcomes.

`ct-validation` tests whether a set of gene-indication pairs is enriched for clinical success. It computes risk ratios and odds ratios with confidence intervals across clinical phase transitions and supports semantic disease matching through ontology-based similarity.

> **Paper:** Kostiuk K, Igumnov D, Fedichev P, Feizi A. _ct-validation: an open framework for benchmarking gene-indication evidence against clinical trial outcomes._ (2026)

## Installation

Requires Python 3.11+.

```bash
pip install ct-validation
```

Optional extras:

```bash
pip install ct-validation[plot]  # forest plot visualization
pip install ct-validation[mcp]   # MCP server for agent workflows
pip install ct-validation[parse] # data source parsers
pip install ct-validation[fetch] # ChEMBL fetching script dependencies
```

To run `notebooks/benchmark.py`, which reproduces the paper, install its pinned environment
instead — Jupyter and the plotting theme on top of an editable install of this package:

```bash
pip install -r notebooks/requirements.txt
```

The full run takes about 10 minutes and peaks at 8.5 GB resident memory, so 16 GB of RAM is
recommended.

## Quick start

### Python API

```python
import ct_validation as ctv

results = ctv.validate(
    clinical_trials="data/clinical_trials/aggregated/gene_indication_max_phase.parquet",
    targets="data/genetic_evidence/aggregated/genetic_evidence.parquet",
    similarity_lookup="data/mappings/efo_similarity_lookup_0.5.parquet",
)
print(results)
#   phase_label  n_yes   n_no  rr  rr_ci_lower  rr_ci_upper  p_value  ...
```

Export annotated trials and/or matched-pair audit rows (returned in order):

```python
enrichment, trials = ctv.validate(..., return_trials=True)
enrichment, matched = ctv.validate(..., return_matched_pairs=True)
enrichment, trials, matched = ctv.validate(..., return_trials=True, return_matched_pairs=True)
```

Or build the match-audit table directly (columns: `gene`, `ct_efo_id`, `ge_efo_id`, `similarity`, `match_type`):

```python
matched = ctv.create_matched_pairs_df(
    genetic_evidence=..., clinical_trials=..., similarity_pairs=...,
    similarity_threshold=0.8,
)
```

Batch mode — compare multiple evidence sources at once:

```python
results = ctv.validate(
    clinical_trials="data/clinical_trials/aggregated/gene_indication_max_phase.parquet",
    targets=[
        "data/genetic_evidence/gwas_catalog.parquet",
        "data/genetic_evidence/clinvar.parquet",
        "data/genetic_evidence/genebass.parquet",
    ],
    similarity_lookup="data/mappings/efo_similarity_lookup_0.5.parquet",
)
# returns a list of DataFrames, one per evidence source
```

Prioritized mode — test whether a novel source adds value over an established baseline:

```python
results = ctv.validate(
    clinical_trials="data/clinical_trials/aggregated/gene_indication_max_phase.parquet",
    targets="data/genetic_evidence/novel_score.parquet",
    baseline_evidence="data/genetic_evidence/established_genetics.parquet",
    similarity_lookup="data/mappings/efo_similarity_lookup_0.5.parquet",
)
# pairs supported only by baseline are excluded
```

Expand a disease set using semantic similarity:

```python
expanded = ctv.get_expanded_disease_set(
    efo_ids={"EFO:0000270", "EFO:0000384"},
    similarity_pairs="data/mappings/efo_similarity_lookup_0.5.parquet",
    similarity_threshold=0.8,
)
```

### CLI

```bash
# With config file
ct-validation --config configs/default.yaml

# With explicit arguments
ct-validation \
    --clinical-trials ct.parquet \
    --targets evidence.parquet \
    --similarity-lookup similarity.parquet \
    -o results/

# Batch mode (multiple evidence sources)
ct-validation \
    --clinical-trials ct.parquet \
    --targets gwas.parquet --targets clinvar.parquet --targets genebass.parquet \
    -o results/

# Save annotated trials and/or matched-pair audit rows
ct-validation --config configs/default.yaml \
    --save-trials --save-matched-pairs -o results/
```

### MCP server

```bash
ct-validation-mcp
```

Exposes two tools for agent-based workflows:

- `ct_validate` — compute phase-transition enrichment (includes `p_value`)
- `expand_disease_set` — expand EFO IDs via semantic similarity

## Input schemas

| Input               | Columns                                                | Description                                                                                 |
| ------------------- | ------------------------------------------------------ | ------------------------------------------------------------------------------------------- |
| `clinical_trials`   | `gene`, `efo_id`, `max_phase`, `is_ongoing` (optional) | Target-indication pairs with highest phase reached, and whether that phase is still running |
| `targets`           | `gene`, `efo_id`                                       | Gene-indication pairs with supporting evidence                                              |
| `similarity_lookup` | `efo_id_1`, `efo_id_2`, `similarity`                   | Pairwise EFO similarity, **symmetric with diagonal** (optional)                             |
| `baseline_evidence` | `gene`, `efo_id`                                       | Baseline evidence for prioritized mode (optional)                                           |
| `gene_universe`     | one gene per line (text file)                          | Restrict analysis to these genes (optional)                                                 |

All inputs accept Parquet files or pandas DataFrames (except `gene_universe`, which is a text file or a Python set).

`clinical_trials` must carry **one row per `(gene, efo_id)` pair** — rows are counted individually, so a repeated pair would inflate the denominators. `validate()` raises on duplicates rather than collapsing them, since only the caller knows what `max_phase` and `is_ongoing` should be for a merged row.

## Output schema

| Column                             | Description                                                 |
| ---------------------------------- | ----------------------------------------------------------- |
| `phase_from`, `phase_to`           | Phase transition (e.g. 1→2, 1→4)                            |
| `n_yes`, `n_no`                    | Pairs entering phase, less censored (with/without evidence) |
| `x_yes`, `x_no`                    | Pairs reaching target phase                                 |
| `rate_yes`, `rate_no`              | Progression rates                                           |
| `rr`, `rr_ci_lower`, `rr_ci_upper` | Risk ratio with 95% CI (Katz log method)                    |
| `or`, `or_ci_lower`, `or_ci_upper` | Odds ratio with 95% CI (Woolf logit method)                 |
| `p_value`                          | Two-sided Fisher's exact test p-value                       |

When either comparison group is empty (`n_yes=0` or `n_no=0`), `rr`, `or`, their confidence intervals, and `p_value` are undefined (`NaN`).

### Clustered bootstrap intervals (optional)

Katz and Woolf assume independent rows, but one gene spans many indications — so their intervals run narrow when the correlated unit is the gene rather than the pair. Pass `cluster_col` to add percentile intervals that resample clusters with replacement instead of rows:

```python
results = ctv.validate(..., cluster_col="gene")
```

This adds four columns — `rr_boot_ci_lower`, `rr_boot_ci_upper`, `or_boot_ci_lower`, `or_boot_ci_upper` — and nothing else: resampling only re-estimates spread, so point estimates, counts and closed-form intervals are unchanged. Without `cluster_col` they are absent, so the default schema is the table above.

`p_value` is unaffected: it stays a row-level Fisher's exact test on the observed table, so it does not
carry the clustering the intervals do.

`bootstrap_replicates` (default 10,000) and `bootstrap_seed` (default 0) control the resampling. The percentile bounds carry Monte Carlo error falling as `1/√n`; at the default each bound shifts by a few percent of the interval's width across seeds. Both are ignored without `cluster_col`.

No continuity correction is applied, unlike the closed-form intervals — a percentile bootstrap needs neither `log(ratio)` nor its standard error. A replicate that empties a cell keeps its ratio, so **a bound may come back `0` or infinite**; relevant if you log-scale these columns.

Two clusters per arm is the floor for a defined interval, not a threshold at which one is
trustworthy — a percentile cluster bootstrap undercovers badly on a few dozen clusters, so read bounds from a small
cluster count as indicative.

Bounds are `NaN` when an arm holds fewer than two clusters, when the observed table has a zero cell, or when over 25% of replicates leave a ratio undefined (`0/0`) — past 5% you get a `BootstrapReplicateLossWarning`. The closed-form columns stay populated unless an arm is empty outright (`n_yes=0` or `n_no=0`).

### Matched-pairs export (optional)

When `return_matched_pairs=True` (API) or `--save-matched-pairs` (CLI) is set, a separate audit table is returned/saved with:

| Column       | Description                                              |
| ------------ | -------------------------------------------------------- |
| `gene`       | Gene symbol                                              |
| `ct_efo_id`  | Disease on the clinical trial row                        |
| `ge_efo_id`  | Supporting genetic-evidence disease                      |
| `similarity` | Match score (`1.0` for exact matches)                    |
| `match_type` | `exact` when `ct_efo_id == ge_efo_id`, else `similarity` |

## Enrichment logic

For each phase transition, target-indication pairs that reached at least the starting phase are divided into supported and unsupported groups. The risk ratio is:

```
RR = (x_yes / n_yes) / (x_no / n_no)
```

A risk ratio greater than one indicates that genetically supported pairs are more likely to progress.

When a similarity lookup is provided, a pair (gene, disease) is considered supported if there exists evidence (gene, disease') with similarity above the threshold (default 0.8). Matching probes a single orientation of each pair, so the lookup **must be symmetric** — every pair `(a, b)` stored as `(b, a)` too, plus a diagonal (self-similarity=1.0). Validate one with `check_similarity_lookup(sim)` (or pass `check_similarity=True` to `validate()`); it raises rather than silently dropping half of each disease neighbourhood.

### Censoring undetermined outcomes

When `clinical_trials` carries an `is_ongoing` column, a pair whose highest phase is still running is censored from any transition it has not yet completed — dropped from both `n` and `x` where `is_ongoing` and `max_phase < phase_to`. Such a pair has not failed the transition; its outcome is not known yet. Censoring is per transition, so a pair that already reached `phase_to` still counts as a success. Without the column every pair counts as concluded, nothing is censored, and `validate()` warns. A missing value in the column is read as not-ongoing, matching that default; non-boolean values raise.

`is_ongoing` is the caller's to define. The parsers in `scripts/` set it where nothing concluded at the pair's highest phase **and** that phase is below approval — reaching approval settles the outcome however many post-marketing trials are still recruiting.

### Prioritized mode

When `baseline_evidence` is provided, pairs supported _only_ by the baseline are excluded. This tests whether a novel evidence source adds predictive value beyond an established benchmark.

## Visualization

```python
import ct_validation as ctv

results = ctv.validate(...)
ctv.forest_plot(results, metric="rr", title="Phase I → Approved")
```

## Data source parsers

The `scripts/` directory contains reproducible parsers for public databases:

**Genetic evidence** (`scripts/parse/genetic_evidence/`):

- GWAS Catalog — genome-wide significant associations (p < 5e-8)
- ClinVar — pathogenic/likely pathogenic variants
- Open Targets — germline genetic evidence streams (score ≥ 0.5); somatic sources
  (`eva_somatic`, `intogen`) are opt-in via `--include-somatic`
- Genebass — exome-wide associations (p < 1e-7)

**Clinical trials** (`scripts/parse/clinical_trials/`):

- ChEMBL — gene-drug links (pChEMBL ≥ 7.0); its indications are parsed but excluded
  from the aggregate, since they carry no trial status to censor on
- Open Targets — known drug and indication data
- STITCH — high-confidence activation/inhibition links
- DGIdb — drug-gene interactions
- TrialPanorama — interventional studies

**Ontology** (`scripts/r/`):

- EFO semantic similarity matrix (Lin + Resnik information content)

**Benchmark comparator** (`scripts/parse/minikel.py`):

- Minikel et al. (2024) — Citeline programmes with per-transition censoring, and their
  genetic associations, for the published-vs-ours comparison

See [DATA_SOURCES.md](DATA_SOURCES.md) for download links, versions, and fetching instructions.

Configure paths in `configs/parsing.yaml` and run:

```bash
python scripts/parse/gene_universe.py
python scripts/parse/minikel.py
python scripts/parse/run_parsing.py
```

## Configuration

See `configs/default.yaml` for validation settings and `configs/parsing.yaml` for data source paths. All config values can be overridden via CLI arguments.

Output options in `configs/default.yaml`:

```yaml
output:
    dir: "results/"
    save_trials: false
    save_matched_pairs: false
```

## Development

CI runs on push/PR to `main` (Python 3.11–3.13, `ruff`, `pytest`):

```bash
uv sync --extra mcp --extra parse
uv run ruff check src tests scripts
uv run pytest tests/ -v
```

## Availability

| Resource                      | Location                                                           |
| ----------------------------- | ------------------------------------------------------------------ |
| Source code                   | https://github.com/gero-science/ct-validation                      |
| Archived code (all versions)  | [10.5281/zenodo.21840865](https://doi.org/10.5281/zenodo.21840865) |
| Benchmark data (all versions) | [10.5281/zenodo.21839216](https://doi.org/10.5281/zenodo.21839216) |
| Package                       | [PyPI](https://pypi.org/project/ct-validation/)                    |

The data deposit carries every input `notebooks/benchmark.py` reads. Extract it at the repository
root and the notebook runs without fetching anything:

```bash
tar xzf ct-validation-data.tar.gz
pip install -r notebooks/requirements.txt
cd notebooks && python benchmark.py
```

See [DATA_SOURCES.md](DATA_SOURCES.md) for what the archive contains, per-source licences, and how
to regenerate anything it omits.

## License

MIT
