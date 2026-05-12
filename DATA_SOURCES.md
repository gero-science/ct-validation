# Data sources

Parsed outputs needed to reproduce the benchmark are included in the data release. This document lists the raw upstream databases used by the parsing scripts.

## Genetic evidence

| Source | Version | Files | Download | Acquisition |
|---|---|---|---|---|
| GWAS Catalog | Release 2026-04-07, Header version v1.0.2 | `gwas-catalog-associations_ontology-annotated-full.zip` | https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/04/07/ | `fetch.sh` |
| ClinVar | 2026-04-04 | `clinvar.vcf.gz` | https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/ | `fetch.sh` |
| Open Targets Platform | 25.12 | `association_by_datasource_direct/`, `target/` | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/ | `fetch.sh` |
| OMIM | 2022-06-09 | `genemap2.txt` | https://omim.org/downloads | Manual (requires OMIM account) |
| Genebass | 500k | `results.mt` | `gs://ukbb-exome-public/500k/results/results.mt` | Manual (requires GCP account); preprocess with `scripts/fetch/genebass_preprocess.py` (requires Hail + Java 11) |

## Clinical trials

| Source | Version | Files | Download | Acquisition |
|---|---|---|---|---|
| ChEMBL | 36 | `chembl_36.db` (SQLite) | https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_36/ | `fetch.sh` (via `chembl_fetch.sh`) |
| Open Targets Platform | 25.12 | `known_drug/` | https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/ | `fetch.sh` |
| STITCH / STRING | v5.0 / v12.0 | `9606.actions.v5.0.tsv.gz`, `chemical_sources.v5.0.tsv.gz`, `9606.protein.aliases.v12.0.txt.gz` (STRING, for gene symbol mapping) | http://stitch.embl.de/cgi/download.pl, https://string-db.org/cgi/download | `fetch.sh` |
| DGIdb | 2024-Dec | `interactions.tsv`, `drugs.tsv` | https://dgidb.org/downloads | `fetch.sh` |
| TrialPanorama | 2025-08-22 | Full dataset | https://huggingface.co/datasets/TrialPanorama/TrialPanorama-database | `fetch.sh` (via `hf download`) |

## Ontology and mappings

| Resource | Version | Files | Download | Acquisition |
|---|---|---|---|---|
| EFO | v3.84.0 | `efo.obo` | https://github.com/EBISPOT/efo/releases | `fetch.sh` |
| OxO cross-references | — | OMIM-to-EFO, MeSH-to-EFO, MONDO/Orphanet/HP-to-EFO | — | Included in data release |
| Pan-UK Biobank phenotype manifest | — | `ukb_phenotype_manifest.tsv` | [Google Sheets](https://docs.google.com/spreadsheets/d/1AeeADtT0U1AukliiNyiVzVRdLYPkTbruQSk38DeutU8) | Included in data release |
| Protein-coding genes | Ensembl 100, GRCh38 | `protein_coding_genes.txt` | — | Included in data release |

## Minikel et al. (2024) benchmark data

The 2x2 benchmark comparison uses genetic evidence and Citeline-derived clinical trial data published by [Minikel et al. (2024)](https://github.com/ericminikel/genetic_support/). The original data uses MeSH disease identifiers; the EFO-remapped versions included in `data/minikel/` were produced by joining with `data/mappings/mesh_to_efo.tsv` (OxO cross-references).

## Parsing

After placing source files according to `configs/parsing.yaml`:

```bash
pip install ct-validation[parse]
python scripts/parse/run_parsing.py
```
