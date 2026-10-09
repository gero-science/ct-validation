#!/bin/bash
#
# Fetch raw data sources for ct_validation.
#
# Usage:
#   ./fetch.sh           # all steps
#   ./fetch.sh 3         # only step 3
#   ./fetch.sh 3-        # step 3 to end
#   ./fetch.sh 3-5       # steps 3..5
#   ./fetch.sh 3,5,7     # specific steps
#   EFO_VERSION=3.90.0 ./fetch.sh 10  # another EFO release (default 3.84.0, the paper's)
#   OT_RELEASE=26.09 ./fetch.sh 3,5   # another Open Targets release (default 25.12, the paper's)

set -eo pipefail

TOTAL=10

declare -A STEP_NAMES=(
    [1]="GWAS Catalog"
    [2]="ClinVar"
    [3]="Open Targets (genetic evidence)"
    [4]="ChEMBL"
    [5]="Open Targets (drugs and indications)"
    [6]="STITCH actions"
    [7]="STITCH chemical sources + STRING protein aliases"
    [8]="DGIdb"
    [9]="TrialPanorama"
    [10]="EFO ontology"
)

CURRENT_STEP=""
on_error() {
    local exit_code=$?
    echo "" >&2
    echo "✗ FAILED at step ${CURRENT_STEP} (exit ${exit_code})" >&2
    exit "$exit_code"
}
trap on_error ERR

WGET="wget -q --show-progress --backups=0"
BASE_DIR="${BASE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
EFO_VERSION="${EFO_VERSION:-3.84.0}"
OT_RELEASE="${OT_RELEASE:-25.12}"

fetch() {
    local url="$1"
    if ! $WGET "$url"; then
        echo "ERROR: failed to fetch ${url}" >&2
        return 1
    fi
}

# Enter the Open Targets release directory, `OT_RELEASE` (YY.MM or `latest`), and keep the
# release's croissant.json there: its `version` records which release `latest` was.
enter_ot_release() {
    local d="${BASE_DIR}/data/sources/opentargets/${OT_RELEASE}"
    mkdir -p "$d" && cd "$d"
    rm -f croissant.json
    fetch "https://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_RELEASE}/croissant.json"
}

# Mirror one Open Targets Platform dataset into the current directory. Spark-written
# datasets carry _SUCCESS; the ones written outside Spark from 26.03 are a single parquet.
fetch_ot() {
    rm -rf "$1"
    $WGET --recursive --no-parent --no-host-directories --cut-dirs 6 --reject "index.html*" \
        "ftp://ftp.ebi.ac.uk/pub/databases/opentargets/platform/${OT_RELEASE}/output/$1/" || true
    [ -f "$1/_SUCCESS" ] || compgen -G "$1/*.parquet" > /dev/null \
        || { echo "ERROR: OT $1 incomplete" >&2; return 1; }
}

step_1() {
    local d="${BASE_DIR}/data/sources/gwas_catalog"
    mkdir -p "$d" && cd "$d"
    rm -f gwas-catalog-associations*.zip* gwas-catalog-associations.tsv*
    fetch https://ftp.ebi.ac.uk/pub/databases/gwas/releases/2026/04/07/gwas-catalog-associations_ontology-annotated-full.zip
    unzip -qo gwas-catalog-associations_ontology-annotated-full.zip
    mv gwas-catalog-download-associations-alt-full.tsv gwas-catalog-associations.tsv
    gzip -f gwas-catalog-associations.tsv
}

step_2() {
    local d="${BASE_DIR}/data/sources/clinvar"
    mkdir -p "$d" && cd "$d"
    rm -f clinvar*.vcf.gz
    fetch https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2026/clinvar_20260404.vcf.gz
    mv clinvar_20260404.vcf.gz clinvar.vcf.gz
}

step_3() {
    enter_ot_release
    fetch_ot association_by_datasource_direct
    # Ontology-propagated table; feeds only the propagated arm in notebooks/benchmark.py §9.
    fetch_ot association_by_datasource_indirect
    # association_overall_direct is not parsed; notebooks/benchmark.py §6a tabulates it as the
    # all-evidence reference row in ot_subsource_comparison.csv
    fetch_ot association_overall_direct
    fetch_ot target
}

step_4() {
    cd "${BASE_DIR}"
    rm -f data/sources/chembl/chembl_36.db data/sources/chembl/chembl_36_sqlite.tar.gz
    rm -rf data/sources/chembl/chembl_36
    bash scripts/fetch/chembl_fetch.sh chembl_36
}

step_5() {
    # The parser reads clinical_report as Open Targets did from 26.06; 26.03 aggregated it
    # under another rule (see DATA_SOURCES.md)
    [[ "$OT_RELEASE" != "26.03" ]] \
        || { echo "ERROR: OT 26.03 clinical data are not supported; use 25.12 or 26.06+" >&2; return 1; }
    enter_ot_release
    # YY.MM releases compare as strings, and `latest` sorts after all of them
    if [[ "$OT_RELEASE" < "26.03" ]]; then
        fetch_ot known_drug
    else
        # known_drug is gone from 26.03: indications come from clinical_report, targets from
        # drug_mechanism_of_action (symbols from target), names from drug_molecule.
        fetch_ot clinical_report
        fetch_ot drug_mechanism_of_action
        fetch_ot drug_molecule
        [ -f target/_SUCCESS ] || fetch_ot target
    fi
}

step_6() {
    local d="${BASE_DIR}/data/sources/stitch"
    mkdir -p "$d" && cd "$d"
    rm -f actions.tsv.gz 9606.actions.v5.0.tsv.gz
    fetch http://stitch.embl.de/download/actions.v5.0/9606.actions.v5.0.tsv.gz
    mv 9606.actions.v5.0.tsv.gz actions.tsv.gz
}

step_7() {
    local d="${BASE_DIR}/data/sources/stitch"
    mkdir -p "$d" && cd "$d"
    rm -f chemical_sources.tsv.gz chemical.sources.v5.0.tsv.gz \
          string_protein_aliases.tsv.gz 9606.protein.aliases.v12.0.txt.gz
    fetch http://stitch.embl.de/download/chemical.sources.v5.0.tsv.gz
    mv chemical.sources.v5.0.tsv.gz chemical_sources.tsv.gz
    fetch https://stringdb-static.org/download/protein.aliases.v12.0/9606.protein.aliases.v12.0.txt.gz
    mv 9606.protein.aliases.v12.0.txt.gz string_protein_aliases.tsv.gz
}

step_8() {
    local d="${BASE_DIR}/data/sources/dgidb"
    mkdir -p "$d" && cd "$d"
    rm -f interactions.tsv* drugs.tsv*
    fetch https://dgidb.org/data/2024-Dec/interactions.tsv
    gzip -f interactions.tsv
    fetch https://dgidb.org/data/2024-Dec/drugs.tsv
    gzip -f drugs.tsv
}

step_9() {
    local d="${BASE_DIR}/data/sources/trialpanorama"
    mkdir -p "$d" && cd "$d"
    # hf download resumes via content-hash cache; no pre-clean
    hf download TrialPanorama/TrialPanorama-database --local-dir . --repo-type dataset
}

step_10() {
    local d="${BASE_DIR}/data/mappings"
    mkdir -p "$d" && cd "$d"
    rm -f efo.obo "efo_v${EFO_VERSION}.obo"
    fetch "https://github.com/EBISPOT/efo/releases/download/v${EFO_VERSION}/efo.obo"
    mv efo.obo "efo_v${EFO_VERSION}.obo"
}

parse_steps() {
    local arg="${1:-}"
    if [[ -z "$arg" ]]; then
        seq 1 "$TOTAL"
    elif [[ "$arg" =~ ^([0-9]+)-([0-9]+)$ ]]; then
        seq "${BASH_REMATCH[1]}" "${BASH_REMATCH[2]}"
    elif [[ "$arg" =~ ^([0-9]+)-$ ]]; then
        seq "${BASH_REMATCH[1]}" "$TOTAL"
    elif [[ "$arg" =~ ^[0-9]+(,[0-9]+)*$ ]]; then
        echo "$arg" | tr ',' '\n'
    else
        echo "Invalid step spec: '$arg' (expected: N | N- | N-M | N,M,...)" >&2
        return 1
    fi
}

main() {
    cd "${BASE_DIR}"
    local steps
    steps=$(parse_steps "${1:-}") || exit 1

    if [[ -z "${1:-}" ]]; then
        echo "[--/${TOTAL}] Genebass: not fetched — the preprocessed associations ship in the data release."
        echo "            Regenerating them from gs://ukbb-exome-public/500k/results/results.mt needs GCP and Hail."
    fi

    for n in $steps; do
        if [[ -z "${STEP_NAMES[$n]:-}" ]]; then
            echo "Unknown step: $n" >&2
            exit 1
        fi
        CURRENT_STEP="${n}/${TOTAL} ${STEP_NAMES[$n]}"
        echo "[${n}/${TOTAL}] ${STEP_NAMES[$n]}..."
        "step_${n}"
    done

    echo "Done."
    if [[ -z "${1:-}" ]]; then
        echo "Mappings (OxO, UKB manifest) and the Genebass associations ship in the data release:"
        echo "  https://doi.org/10.5281/zenodo.21839216"
        echo "Gene universe: python scripts/parse/gene_universe.py"
    fi
}

main "$@"
