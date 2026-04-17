#!/bin/bash

set -e

TOTAL=10
step() { echo "[${1}/${TOTAL}] ${2}..."; }
WGET="wget -q --show-progress --backups=0"

BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
cd "${BASE_DIR}"

mkdir -p ./data/sources/gwas_catalog
mkdir -p ./data/sources/clinvar
mkdir -p ./data/sources/opentargets/25.12
mkdir -p ./data/sources/omim
mkdir -p ./data/sources/genebass
mkdir -p ./data/sources/stitch
mkdir -p ./data/sources/chembl
mkdir -p ./data/sources/dgidb
mkdir -p ./data/sources/trialpanorama
mkdir -p ./data/mappings

echo "[--/${TOTAL}] OMIM: skipped (manual — download genemap2.txt from omim.org/downloads)"
echo "[--/${TOTAL}] Genebass: skipped (manual — requires GCP account for gs://ukbb-exome-public/500k/results/results.mt)"

step 1 "GWAS Catalog"
cd "${BASE_DIR}/data/sources/gwas_catalog"
$WGET https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/04/07/gwas-catalog-associations_ontology-annotated-full.zip
unzip -qo gwas-catalog-associations_ontology-annotated-full.zip
mv gwas-catalog-download-associations-alt-full.tsv gwas-catalog-associations.tsv
gzip -f gwas-catalog-associations.tsv

step 2 "ClinVar"
cd "${BASE_DIR}/data/sources/clinvar"
$WGET https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/clinvar_20260404.vcf.gz
mv clinvar_20260404.vcf.gz clinvar.vcf.gz

step 3 "Open Targets (genetic evidence)"
cd "${BASE_DIR}/data/sources/opentargets/25.12"
$WGET --recursive --no-parent --no-host-directories --cut-dirs 6 --reject "index.html*" ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/association_by_datasource_direct/ || true
$WGET --recursive --no-parent --no-host-directories --cut-dirs 6 --reject "index.html*" ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/target/ || true

step 4 "ChEMBL"
cd "${BASE_DIR}"
chmod +x scripts/fetch/chembl_fetch.sh
scripts/fetch/chembl_fetch.sh chembl_36

step 5 "Open Targets (known drugs)"
cd "${BASE_DIR}/data/sources/opentargets/25.12"
$WGET --recursive --no-parent --no-host-directories --cut-dirs 6 --reject "index.html*" ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/25.12/output/known_drug/ || true

step 6 "STITCH actions"
cd "${BASE_DIR}/data/sources/stitch"
$WGET http://stitch.embl.de/download/actions.v5.0/9606.actions.v5.0.tsv.gz
mv 9606.actions.v5.0.tsv.gz actions.tsv.gz

step 7 "STITCH chemical sources + STRING protein aliases"
$WGET http://stitch.embl.de/download/chemical.sources.v5.0.tsv.gz
mv chemical.sources.v5.0.tsv.gz chemical_sources.tsv.gz
$WGET https://stringdb-static.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz
mv 9606.protein.aliases.v12.0.txt.gz string_protein_aliases.tsv.gz

step 8 "DGIdb"
cd "${BASE_DIR}/data/sources/dgidb"
$WGET https://dgidb.org/data/2024-Dec/interactions.tsv
gzip -f interactions.tsv
$WGET https://dgidb.org/data/2024-Dec/drugs.tsv
gzip -f drugs.tsv

step 9 "TrialPanorama"
cd "${BASE_DIR}/data/sources/trialpanorama"
hf download TrialPanorama/TrialPanorama-database --local-dir . --repo-type dataset

step 10 "EFO ontology"  # bonus step — not counted in skips
cd "${BASE_DIR}/data/mappings"
$WGET https://github.com/EBISPOT/efo/releases/download/v3.84.0/efo.obo
mv efo.obo efo_v3.84.0.obo

echo "Done. Mappings (OxO, UKB manifest, gene universe) must be placed manually from the data release."
