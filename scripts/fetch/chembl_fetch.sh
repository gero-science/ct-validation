#!/usr/bin/env bash
# Download ChEMBL SQLite database from EBI FTP.
#
# Usage:
#   bash scripts/fetch/chembl_fetch.sh              # defaults to chembl_36
#   bash scripts/fetch/chembl_fetch.sh chembl_35    # specific version
set -euo pipefail

CHEMBL_VERSION="${1:-chembl_36}"
BASE_URL="https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases"
OUTPUT_DIR="data/sources/chembl"

mkdir -p "$OUTPUT_DIR"

TARBALL="${CHEMBL_VERSION}_sqlite.tar.gz"
URL="${BASE_URL}/${CHEMBL_VERSION}/${TARBALL}"

echo "Downloading ${CHEMBL_VERSION} SQLite database..."
curl -L --progress-bar -o "${OUTPUT_DIR}/${TARBALL}" "$URL"

echo "Extracting..."
tar -xzf "${OUTPUT_DIR}/${TARBALL}" -C "$OUTPUT_DIR"

# Move DB to a cleaner path and remove nested directories
# Tarball extracts as: chembl_36/chembl_36_sqlite/chembl_36.db
mv "${OUTPUT_DIR}/${CHEMBL_VERSION}/${CHEMBL_VERSION}_sqlite/${CHEMBL_VERSION}.db" "${OUTPUT_DIR}/${CHEMBL_VERSION}.db"
rm -rf "${OUTPUT_DIR}/${CHEMBL_VERSION}" "${OUTPUT_DIR}/${TARBALL}"

DB_PATH="${OUTPUT_DIR}/${CHEMBL_VERSION}.db"
echo "Done: ${DB_PATH}"
echo "Size: $(du -h "${DB_PATH}" | cut -f1)"
