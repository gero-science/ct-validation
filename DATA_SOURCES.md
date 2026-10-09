# Data sources

This document lists the raw upstream databases the parsing scripts consume, and what the
data release ships.

## Data release

The release archive ([10.5281/zenodo.21839216](https://doi.org/10.5281/zenodo.21839216)) carries
everything `notebooks/benchmark.py` reads, so the benchmark runs without re-fetching anything:

| Path                              | Contents                                                                                      |
| --------------------------------- | --------------------------------------------------------------------------------------------- |
| `data/genetic_evidence/`          | the four parsed sources, `opentargets_with_somatic.parquet`, and both aggregates              |
| `data/clinical_trials/`           | the five parsed sources, `aggregated/`, and `aggregated_unknown_ongoing/`                     |
| `data/minikel/`                   | the three processed parquets                                                                  |
| `data/mappings/`                  | similarity lookup at 0.5, gene universe, `efo.obo`, `d2023.bin`, and the crosswalks           |
| `data/sources/opentargets/25.12/` | `target/`, `known_drug/`, and the three association tables                                    |
| `data/sources/dgidb/`             | `drugs.tsv.gz`                                                                                |
| `data/sources/genebass/`          | `gene_associations_p1e-7.parquet`, the Hail-preprocessed intermediate                         |

The Open Targets and DGIdb entries under `sources/` are verbatim upstream releases rather than
our derivatives; they ship because the notebook opens them directly — `known_drug/` for the Open
Targets trial arm, the association tables for the per-subsource breakdown and the prior-art
comparison. The Genebass file is ours: the output of `genebass_preprocess.py`, included so the
parsing step can be re-run without a GCP account or Hail.

Deliberately excluded:

- **`efo_full_similarity_matrix.parquet`** (2.5 GB). `scripts/r/create_efo_similarity_matrix.R`
  writes it and derives the 0.5 lookup from it; nothing reads it afterwards, and every analysis
  here thresholds at 0.5 or above. Regenerate it from that script if you need lower.
- **Raw Citeline files** (`pp.tsv`, `assoc.tsv.gz`, `indic.tsv`). Minikel et al. redistribute
  these themselves; take them from [their repository](https://github.com/ericminikel/genetic_support/).
- **Everything else under `data/sources/`.** Versioned upstream releases; `fetch.sh` retrieves them.

## Open Targets and EFO releases

By default `fetch.sh` downloads Open Targets 26.09 and EFO 3.93, the release Open Targets built
26.09 with, and the configs read them: `parsing.yaml` the four 26.09 clinical datasets and
`efo_v3.93.0.obo`, and `default.yaml` the EFO 3.93 similarity lookup, which
`create_efo_similarity_matrix.R` builds by default (see [Similarity lookup](#similarity-lookup)).
`OT_RELEASE` (YY.MM) and `EFO_VERSION` select others; `EFO_VERSION=ot` reads the EFO release
`OT_RELEASE` was built with from its `manifest.json`, which records it from 25.12. Releases from
25.12 come from Open Targets' public S3 bucket, which serves the same files several times faster
than the FTP but keeps recent releases only; earlier ones come from the FTP.

A release's disease ids follow its own EFO, so the lookup must be built from that EFO: ids EFO
has replaced since find no match. The paper used Open Targets 25.12, EFO 3.84 and no remapping.
To parse it again, fetch them:

```bash
OT_RELEASE=25.12 EFO_VERSION=3.84.0 ./scripts/fetch/fetch.sh 3,5,10
```

then, in `configs/parsing.yaml`, point the Open Targets paths at `25.12`, use
`opentargets_known_drug` (commented there) in place of the four clinical datasets, and set
`efo_ontology` to `efo_v3.84.0.obo` and `efo_remap_obsolete` to false; in `configs/default.yaml`,
use the data release's `efo_similarity_lookup_0.5.parquet`. Every other source is pinned to the
version listed below. `notebooks/benchmark.py` reads the data release directly, so it reproduces
the paper without fetching anything.

## Genetic evidence

| Source                                      | Version                           | Files                                                                         | Download                                                               | Acquisition                                                                                                                                                                             | Licence                                                                                                                          |
| ------------------------------------------- | --------------------------------- | ----------------------------------------------------------------------------- | ---------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- |
| GWAS Catalog                                | Release 2026-04-07, header v1.0.2 | `gwas-catalog-associations_ontology-annotated-full.zip`                       | https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/04/07/          | `fetch.sh`                                                                                                                                                                              | None named. [EMBL-EBI terms](https://www.ebi.ac.uk/about/terms-of-use) add no restrictions beyond the data owners'               |
| ClinVar                                     | 2026-04-04                        | `clinvar.vcf.gz`                                                              | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/  | `fetch.sh`                                                                                                                                                                              | None named. US public domain; NCBI places no restrictions on distribution, though submitters may hold rights                     |
| Open Targets Platform                       | 26.09 (paper: 25.12)              | `association_by_datasource_direct/`, `association_overall_direct/`, `target/` | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/ | `fetch.sh` (S3 from 25.12)                                                                                                                                                              | CC0 1.0                                                                                                                          |
| Open Targets Platform (ontology-propagated) | 26.09 (paper: 25.12)              | `association_by_datasource_indirect/`                                         | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/ | `fetch.sh` (S3 from 25.12)                                                                                                                                                              | CC0 1.0                                                                                                                          |
| Genebass                                    | 500k                              | `results.mt`                                                                  | `gs://ukbb-exome-public/500k/results/results.mt`                       | Manual (GCP account); preprocess with `scripts/fetch/genebass_preprocess.py` (Hail + Java 11). The preprocessed output is in the data release, so neither is needed to re-run the parse | **None stated.** Public bulk downloads of UK Biobank-derived summary statistics; GCP/Hail is an access burden, not a restriction |

Only Open Targets grants a licence. NCBI and EMBL-EBI state that _they_ add no restrictions of
their own, which is not the same as a grant, and Genebass names nothing — so a data-availability
statement should not describe any of these three as CC0.

`association_overall_direct/` is not parsed; it supplies the all-evidence reference row in
`ot_subsource_comparison.csv`.

From 26.03 the association tables name the datasource in `aggregationValue` (under
`aggregationType` `datasourceId`) and the score in `associationScore`; `opentargets.py` reads
either schema to the same output.

Somatic datasources (`eva_somatic`, `intogen`) are excluded by default, matching Minikel et al.'s
germline restriction, and restored by `opentargets.py --include-somatic` for the scope
sensitivity in §3 of the notebook, which writes `opentargets_with_somatic.parquet` alongside the
germline `opentargets.parquet` rather than replacing it. The other three sources need no such filter: ClinVar excludes
somatic conditions by pattern, and GWAS Catalog and Genebass are germline by construction.

The pipeline default is the direct table — evidence attached to a disease term itself, not
inherited from its subclasses. The indirect table rolls the same columns up the disease ontology
(4.7M rows direct against 14.2M indirect) and is read only by the prior-art comparison in
`benchmark.py` §9. Propagation moves evidence toward the ontology root, so most pairs that arm
adds match on near-root terms such as `EFO:0000651` (phenotype), and are enriched for approval
(83% against a 60% base rate) — it measures what naive propagation does to this benchmark, not a
better disease match.

## Clinical trials

| Source                | Version      | Files                                                                                           | Download                                                                  | Acquisition                        |
| --------------------- | ------------ | ----------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------- | ---------------------------------- |
| ChEMBL                | 36           | `chembl_36.db` (SQLite)                                                                         | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/   | `fetch.sh` (via `chembl_fetch.sh`) |
| Open Targets Platform | 26.09        | `clinical_report/`, `drug_mechanism_of_action/`, `drug_molecule/`; paper (25.12): `known_drug/` | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/26.09/output/    | `fetch.sh` (S3 from 25.12)         |
| STITCH / STRING       | v5.0 / v12.0 | `9606.actions.v5.0.tsv.gz`, `chemical_sources.v5.0.tsv.gz`, `9606.protein.aliases.v12.0.txt.gz` | http://stitch.embl.de/cgi/download.pl, https://string-db.org/cgi/download | `fetch.sh`                         |
| DGIdb                 | 2024-Dec     | `interactions.tsv`, `drugs.tsv`                                                                 | https://dgidb.org/downloads                                               | `fetch.sh`                         |
| TrialPanorama         | 2025-08-22   | Full dataset                                                                                    | https://huggingface.co/datasets/TrialPanorama/TrialPanorama-database      | `fetch.sh` (via `hf download`)     |

_Open Targets from 26.06._ `known_drug` was dropped from the Platform at 26.03. From 26.06
`opentargets.py` writes the same two outputs from `clinical_report/`,
`drug_mechanism_of_action/`, `drug_molecule/` and `target/`, each its own entry in
`configs/parsing.yaml` in place of `opentargets_known_drug`. It applies the quality controls
Open Targets uses to aggregate reports into its own `clinical_indication`: reports flagged
`INDIRECT_PRIMARY_PURPOSE` or `PHASE_IV_NOT_APPROVED` are dropped, and `UNVALIDATED_INDICATION`
is only a warning. Under those controls alone, its drug–disease pairs reproduce
`clinical_indication` exactly at 26.06 (86,468) and 26.09 (86,451). The parser then keeps only
interventional registry trials (84,102 and 84,273) and drops reports that establish no clinical
phase (`IND`, `PRECLINICAL`, `UNKNOWN`) or whose trial was withdrawn, so it outputs 74,071 and
74,325 pairs.

26.03 is not supported, and `fetch.sh` refuses it at step 5: Open Targets then dropped
`UNVALIDATED_INDICATION` reports from `clinical_indication`, a rule it reversed at 26.06, and its
`drug_molecule` stores synonyms and trade names as plain strings. Its association tables are
read.

### Phase, status and licensing per source

The internal phase scale is 0 = Preclinical … 4 = Approved.

| Source                    | Contributes | Phase handling                                                   | Study-type filter                        | Status                              | Licence                          |
| ------------------------- | ----------- | ---------------------------------------------------------------- | ---------------------------------------- | ----------------------------------- | -------------------------------- |
| Open Targets `known_drug` | both legs   | Numeric; phases 1–4 kept, the 6,562 records at phase 0.5 dropped | None available at 25.12                  | `status`                            | CC0 1.0                          |
| Open Targets from 26.06   | both legs   | `clinicalStage`, standardised upward                             | Interventional only                      | `trialOverallStatus`                | CC0 1.0                          |
| TrialPanorama             | both legs   | String, standardised upward (see below)                          | Interventional only; own conditions only | `recruitment_status`                | Apache 2.0                       |
| ChEMBL                    | gene–drug   | pChEMBL ≥ 7.0; indications excluded from the aggregate           | n/a                                      | None recorded — hence the exclusion | CC BY-SA 3.0                     |
| STITCH                    | gene–drug   | n/a                                                              | n/a                                      | n/a                                 | CC BY 4.0 (see below)            |
| DGIdb                     | gene–drug   | n/a                                                              | n/a                                      | n/a                                 | MIT code, mixed data (see below) |

_TrialPanorama phases._ Counted after the interventional filter: `EARLY_PHASE1` → `PHASE1`
(5,473 studies); combined phases round up to their highest component (`PHASE1/PHASE2` → `PHASE2`,
15,269; `PHASE2/PHASE3` → `PHASE3`, 7,002). `PHASE1/PHASE2/PHASE3` maps to `PHASE3` and
`PHASE1/PHASE4` matches no rule, but neither contributes — all 11,735 and 770 carry no
`trial_type`. Studies with no phase are dropped (197,071 interventional, 1,072,669 overall).

_TrialPanorama filters._ `trial_type == INTERVENTIONAL` keeps 401,217 of 1,332,141 studies.
Conditions are restricted to `condition_mesh_type == mesh-list`, the study's own conditions,
dropping the 6,230,140 of 8,051,812 rows that are TrialPanorama's ancestor closure over the MeSH
tree — which would otherwise enter a trial of one carcinoma as a trial of _Neoplasms_.

_STITCH licence._ CC BY 4.0 covers every file used (`actions`, `chemical_sources`,
`protein_chemical_links`) and the STRING alias file. The `*.detailed` and `*.transfer` files are
CC BY-NC-SA 4.0 and are not used; three SQL dumps need a separate licence.

_DGIdb licence._ Code is MIT, but the aggregated interactions inherit their upstream sources'
terms, some non-commercial. `dgidb.py` does not retain `interaction_source_db_name`, so upstream
sources cannot be separated from the parsed output. The column is in the raw `interactions.tsv`
of the pinned 2024-Dec release that `fetch.sh` retrieves, so re-parse from there if a
licence-filtered subset is needed.
The terms are an inherited ambiguity rather than a stated bar, which is why DGIdb is retained;
a source whose own licence prohibited commercial use would not be.

_Status._ Classified identically across vocabularies by `trial_status.py`: only `Withdrawn` fails
to establish that a phase was reached; `Completed`, `Terminated`, `Unknown status` and NULL all
count as concluded. [CT.gov](https://clinicaltrials.gov/) applies `Unknown` only once a
completion date has passed, so these went dark rather than stalling. `aggregate.py --unknown-ongoing` produces the sensitivity arm
that censors them instead. NULL is concluded in `known_drug`, where only approval records lack a
status. From 26.06 curated resources lack one too so a status-less report concludes only at `APPROVAL` or
`WITHDRAWAL`, or when `phaseFromSource` says discontinued, and is otherwise left open.

**Two asymmetries between the status-bearing sources are known and deliberately not fixed at
25.12.** Both are gone from 26.06, whose reports carry the study type and keep early phase I.

_Interventional filter._ TrialPanorama is restricted to interventional studies; `known_drug`
carries no study-type field at 25.12, so its leg admits records TrialPanorama's filter would have
removed. Fixing this needs a field the release does not contain.

_Early phase._ Open Targets drops phase-0.5 records while TrialPanorama folds `EARLY_PHASE1` into
Phase I, so the same study can be absent from one leg and present as Phase I in the other. Since
`max_phase` is a maximum across sources, this can only raise a pair into the Phase I cohort,
never lower it out of a higher one.

## Ontology and mappings

| Resource                          | Version            | Files                                 | Download                                                                                             | Acquisition                                                                                                                                         |
| --------------------------------- | ------------------ | ------------------------------------- | ---------------------------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| EFO                               | v3.93.0 (paper: v3.84.0) | `efo_v3.93.0.obo`               | https://github.com/EBISPOT/efo/releases                                                              | `fetch.sh` (`EFO_VERSION`)                                                                                                                          |
| MeSH descriptors                  | 2023               | `d2023.bin`                           | https://nlmpubs.nlm.nih.gov/projects/mesh/2023/asciimesh/                                            | Data release; read only by the benchmark notebook (oncology flag, tree C04)                                                                         |
| OxO cross-references              | OxO 1, retrieved on or before 2026-05-08 | MeSH-to-EFO, MONDO/Orphanet/HP-to-EFO | https://www.ebi.ac.uk/spot/oxo                                                                       | Data release, which is the reference copy: a dated snapshot of a live service, so a later query is not guaranteed to reproduce it                    |
| Pan-UK Biobank phenotype manifest | —                  | `ukb_phenotype_manifest.tsv`          | [Google Sheets](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8) | Data release                                                                                                                                        |
| Protein-coding genes              | Open Targets 26.09 | `protein_coding_genes.txt`            | —                                                                                                    | `python scripts/parse/gene_universe.py` (`target.biotype == protein_coding`; 20,080 symbols, 624 of them Ensembl IDs for genes with no HGNC symbol; 20,097 and 634 at 25.12, the paper's) |

### EFO releases and obsolete terms

Each source is coded against the EFO release it was built with, and the similarity lookup
against one release. Matching joins through the lookup, so an `efo_id` the lookup lacks can never
match: `validate()` leaves clinical-trial pairs with one out of both groups, and evidence with
one supports nothing.

With `mappings.efo_remap_obsolete` in `configs/parsing.yaml`, **every
parser rewrites obsolete ids to their live successor** in `mappings.efo_ontology`, following
`replaced_by`, before the source groups by disease; two old ids replaced by one term
therefore collapse into one pair with the right `max_phase` and `is_ongoing`. It is on
by default; the paper's parsed data were built with it off.

Against EFO 3.93, the release Open Targets 26.09 was built with, this leaves no obsolete ids in
any source. What it cannot fix is terms EFO does not import: about 8–9% of TrialPanorama and
Minikel CT records, nearly all Orphanet ids from the MeSH crosswalk. Neither does it make all
evidence matchable: the R script leaves the OBA (measurements) and PR (proteins) branches out of the lookup
by design (`--exclude-prefixes`), in the paper as now, and OBA terms make up 14–24% of GWAS
Catalog, Genebass and Open Targets genetic evidence records.

_ClinVar._ Its third mapping tier queries OLS live, so its hits depend on the date of the query,
not on the configured EFO release; they are cached in `clinvar_ols_cache.csv`, which a re-run
reuses. With remapping on, all three tiers accept only terms live in the configured release or
with a live successor, and the cached hits are remapped like the rest.

### Similarity lookup

`scripts/r/create_efo_similarity_matrix.R` builds the EFO similarity lookup. Its R environment
is pinned in `scripts/r/pixi.toml` (install [pixi](https://pixi.sh) first; on Apple Silicon,
`ontologySimilarity` is compiled from CRAN on the first run). From the repository root:

```bash
pixi run --manifest-path scripts/r/pixi.toml similarity
```

By default it reads `efo_v3.93.0.obo` and writes `efo_similarity_lookup_v3.93.0_0.5.parquet`,
which `configs/default.yaml` reads, and the full matrix next to it; `--efo-path`, `--output` and
`--lookup-output` select another release. EFO 3.93 takes about 22 minutes and 5 GB of RAM on one
core.

## Minikel et al. (2024) benchmark data

The comparison uses genetic evidence and Citeline-derived trial data published by
[Minikel et al. (2024)](https://github.com/ericminikel/genetic_support/). Their data uses MeSH;
`scripts/parse/minikel.py` remaps both legs to EFO via `mesh_to_efo.tsv` and writes to
`data/minikel/`. The raw files are a manual download from their repository.

| Resource                   | Files                                | Download                                                                   |
| -------------------------- | ------------------------------------ | -------------------------------------------------------------------------- |
| Pharmaprojects programmes  | `pp.tsv`                             | https://github.com/ericminikel/genetic_support/tree/main/data              |
| Genetic associations       | `assoc.tsv.gz`                       | same                                                                       |
| Indication genetic insight | `indic.tsv`                          | same                                                                       |
| Published risk ratios      | `minikel_enrichment_results.parquet` | Data release; used only for the sanity check against their reported values |

_Censoring._ `minikel.py` derives `is_ongoing` from the Pharmaprojects active category (`acat`):
a programme with an active category has concluded nothing at its current phase. This reproduces
their censoring exactly — their per-transition `succ_*` flags are NA for precisely the programmes
this rule censors, with zero disagreements over 29,313 programmes.

_Associations._ Drops `intOGen` (somatic, matching the germline restriction above) and keeps Open
Targets Genetics rows only at `l2g_share >= 0.5`, those being locus-to-gene predictions rather
than gene-level calls. Gene-level sources carry no `l2g_share` and are unaffected.

_Denominator._ Indications with `genetic_insight == "none"` are dropped when `indic.tsv` is
configured, matching their restriction to diseases where genetic evidence could exist. Their
published ratios are computed after this: reported cohorts are 13,022 / 7,223 / 2,184 for the
three transitions, against 13,030 / 7,225 / 2,184 with the filter and 14,770 / 8,340 / 2,604
without. Our own pipeline applies no such filter, so the cross-source comparison also carries a
difference in which diseases are eligible.

Their composite I→Launch statistic is _not_ re-censored: its denominator is the phase-I cohort,
so programmes still running at Phase II or III count as non-launches. Per-transition censoring of
I→Approved therefore gives a different quantity from their published 2.63.

## Parsing

After placing source files according to `configs/parsing.yaml`:

```bash
pip install ct-validation[parse]
python scripts/parse/gene_universe.py
python scripts/parse/minikel.py
python scripts/parse/run_parsing.py
```
