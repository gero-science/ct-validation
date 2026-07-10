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

## Quick start

### Python API

```python
import ct_validation as ctv

results = ctv.validate(
    clinical_trials="data/clinical_trials/gene_indication_max_phase.parquet",
    targets="data/genetic_evidence/genetic_evidence.parquet",
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
    clinical_trials="data/clinical_trials/gene_indication_max_phase.parquet",
    targets=[
        "data/genetic_evidence/gwas_catalog.parquet",
        "data/genetic_evidence/clinvar.parquet",
        "data/genetic_evidence/omim.parquet",
    ],
    similarity_lookup="data/mappings/efo_similarity_lookup_0.5.parquet",
)
# returns a list of DataFrames, one per evidence source
```

Prioritized mode — test whether a novel source adds value over an established baseline:

```python
results = ctv.validate(
    clinical_trials="data/clinical_trials/gene_indication_max_phase.parquet",
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
    --targets gwas.parquet --targets clinvar.parquet --targets omim.parquet \
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

| Input               | Columns                              | Description                                        |
| ------------------- | ------------------------------------ | -------------------------------------------------- |
| `clinical_trials`   | `gene`, `efo_id`, `max_phase`        | Target-indication pairs with highest phase reached |
| `targets`           | `gene`, `efo_id`                     | Gene-indication pairs with supporting evidence     |
| `similarity_lookup` | `efo_id_1`, `efo_id_2`, `similarity` | Pairwise EFO similarity, **symmetric with diagonal** (optional) |
| `baseline_evidence` | `gene`, `efo_id`                     | Baseline evidence for prioritized mode (optional)  |
| `gene_universe`     | one gene per line (text file)         | Restrict analysis to these genes (optional)        |

All inputs accept Parquet files or pandas DataFrames (except `gene_universe`, which is a text file or a Python set).

## Output schema

| Column                             | Description                                  |
| ---------------------------------- | -------------------------------------------- |
| `phase_from`, `phase_to`           | Phase transition (e.g. 1→2, 1→4)             |
| `n_yes`, `n_no`                    | Pairs entering phase (with/without evidence) |
| `x_yes`, `x_no`                    | Pairs reaching target phase                  |
| `rate_yes`, `rate_no`              | Progression rates                            |
| `rr`, `rr_ci_lower`, `rr_ci_upper` | Risk ratio with 95% CI (Katz log method)     |
| `or`, `or_ci_lower`, `or_ci_upper` | Odds ratio with 95% CI (Woolf logit method)  |
| `p_value`                          | Two-sided Fisher's exact test p-value        |

When either comparison group is empty (`n_yes=0` or `n_no=0`), `rr`, `or`, their confidence intervals, and `p_value` are undefined (`NaN`).

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

A risk ratio greater than one indicates that genetically supported pairs are more likely to progress. When a similarity lookup is provided, a pair (gene, disease) is considered supported if there exists evidence (gene, disease') with similarity above the threshold (default 0.8). Matching probes a single orientation of each pair, so the lookup **must be symmetric** — every pair `(a, b)` stored as `(b, a)` too, plus a diagonal (self-similarity=1.0). Validate one with `check_similarity_lookup(sim)` (or pass `check_similarity=True` to `validate()`); it raises rather than silently dropping half of each disease neighbourhood.

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

- GWAS Catalog — genome-wide significant associations (p < 1e-8)
- ClinVar — pathogenic/likely pathogenic variants
- OMIM — established molecular basis (mapping code 3)
- Open Targets — genetic evidence streams (score ≥ 0.5)
- Genebass — exome-wide associations (p ≤ 1e-7)

**Clinical trials** (`scripts/parse/clinical_trials/`):

- ChEMBL — gene-drug and drug-indication links (pChEMBL > 7.0)
- Open Targets — known drug and indication data
- STITCH — high-confidence activation/inhibition links
- DGIdb — drug-gene interactions
- TrialPanorama — interventional studies

**Ontology** (`scripts/r/`):

- EFO semantic similarity matrix (Lin + Resnik information content)

See [DATA_SOURCES.md](DATA_SOURCES.md) for download links, versions, and fetching instructions.

Configure paths in `configs/parsing.yaml` and run:

```bash
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
uv sync --extra mcp
uv run ruff check src tests
uv run pytest tests/ -v
```

## License

MIT
