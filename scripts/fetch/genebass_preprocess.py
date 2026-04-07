#!/usr/bin/env python3
"""Preprocess Genebass results.mt to filtered parquet.

Reads the raw Hail MatrixTable from Genebass (500k exome sequencing)
and exports significant gene-phenotype associations as a parquet file.

Requires Hail (pip install hail) and Java 17.

Input:
    results.mt from gs://ukbb-exome-public/500k/results/results.mt

Output:
    gene_associations_p1e-5.parquet with columns:
        gene_id, gene_symbol, annotation, phenocode, description,
        coding_description, category, n_cases, n_controls,
        Pvalue, Pvalue_Burden, Pvalue_SKAT, BETA_Burden, SE_Burden
"""

import argparse
from pathlib import Path

import hail as hl


def preprocess_genebass(
    input_path: str,
    output_path: str,
    pvalue_threshold: float = 1e-5,
    n_cores: int = 16,
) -> None:
    """Filter Genebass MatrixTable by p-value and export to parquet."""
    hl.init(master=f"local[{n_cores}]", log="hail.log", quiet=True)

    print(f"Loading MatrixTable from: {input_path}")
    mt = hl.read_matrix_table(input_path)

    n_genes = mt.count_rows()
    n_phenotypes = mt.count_cols()
    print(f"Genes: {n_genes:,}, Phenotypes: {n_phenotypes:,}")

    print(f"Filtering p-value < {pvalue_threshold}")
    mt_filtered = mt.filter_entries(mt.Pvalue < pvalue_threshold)

    results_table = mt_filtered.entries().key_by()
    results_table = results_table.select(
        "gene_id",
        "gene_symbol",
        "annotation",
        "phenocode",
        "description",
        "coding_description",
        "category",
        "n_cases",
        "n_controls",
        "Pvalue",
        "Pvalue_Burden",
        "Pvalue_SKAT",
        "BETA_Burden",
        "SE_Burden",
    )

    n_results = results_table.count()
    print(f"Significant associations: {n_results:,}")

    print(f"Exporting to: {output_path}")
    results_table.export(output_path)
    print("Done!")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--input",
        type=str,
        default="data/sources/genebass/results.mt",
        help="Path to Genebass results.mt",
    )
    parser.add_argument(
        "--output",
        type=str,
        default="data/sources/genebass/gene_associations_p1e-5.parquet",
        help="Output parquet path",
    )
    parser.add_argument(
        "--pvalue",
        type=float,
        default=1e-5,
        help="P-value threshold (default: 1e-5)",
    )
    parser.add_argument(
        "--cores",
        type=int,
        default=16,
        help="Number of cores for Hail (default: 16)",
    )
    args = parser.parse_args()

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    preprocess_genebass(args.input, args.output, args.pvalue, args.cores)


if __name__ == "__main__":
    main()
