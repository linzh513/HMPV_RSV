#!/usr/bin/env python3

import argparse
import sys
import time
from functools import partial
from multiprocessing import Pool, cpu_count
from pathlib import PurePath

import pandas as pd
from arboreto.algo import _prepare_input
from arboreto.core import (
    EARLY_STOP_WINDOW_LENGTH,
    RF_KWARGS,
    SGBM_KWARGS,
    infer_partial_network,
    target_gene_indices,
    to_tf_matrix,
)
from arboreto.utils import load_tf_names
from pyscenic.cli.utils import load_exp_matrix, suffixes_to_separator
from tqdm import tqdm


def create_argument_parser():
    """Create the command-line argument parser."""
    parser = argparse.ArgumentParser(
        description="Run Arboreto using a multiprocessing pool."
    )

    parser.add_argument(
        "expression_mtx_fname",
        type=str,
        help=(
            "Expression matrix file. Supported formats are CSV and loom. "
            "For CSV files, rows should be cells and columns should be genes. "
            "For loom files, rows should be genes and columns should be cells."
        ),
    )

    parser.add_argument(
        "tfs_fname",
        type=str,
        help="Text file containing transcription factor names, one TF per line.",
    )

    parser.add_argument(
        "-m",
        "--method",
        choices=["genie3", "grnboost2"],
        default="grnboost2",
        help="Network reconstruction method. Default: grnboost2.",
    )

    parser.add_argument(
        "-o",
        "--output",
        type=str,
        default="-",
        help="Output TSV or CSV file. Use '-' to write to standard output.",
    )

    parser.add_argument(
        "--num_workers",
        type=int,
        default=cpu_count(),
        help=f"Number of worker processes. Default: {cpu_count()}.",
    )

    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Optional random seed for regressor initialization.",
    )

    parser.add_argument(
        "--cell_id_attribute",
        type=str,
        default="CellID",
        help="Cell identifier attribute in the loom file.",
    )

    parser.add_argument(
        "--gene_attribute",
        type=str,
        default="Gene",
        help="Gene identifier attribute in the loom file.",
    )

    parser.add_argument(
        "--sparse",
        action="store_true",
        help="Load the expression matrix as a sparse CSC matrix.",
    )

    parser.add_argument(
        "-t",
        "--transpose",
        action="store_true",
        help="Transpose the expression matrix.",
    )

    return parser


def run_infer_partial_network(
    target_gene_index,
    gene_names,
    expression_matrix,
    tf_matrix,
    tf_matrix_gene_names,
    method_params,
    seed,
):
    """Infer regulatory links for one target gene."""
    target_gene_name = gene_names[target_gene_index]
    target_gene_expression = expression_matrix[:, target_gene_index]

    return infer_partial_network(
        regressor_type=method_params[0],
        regressor_kwargs=method_params[1],
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        target_gene_name=target_gene_name,
        target_gene_expression=target_gene_expression,
        include_meta=False,
        early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
        seed=seed,
    )


def get_method_parameters(method):
    """Return the regressor type and parameters for the selected method."""
    if method == "grnboost2":
        return "GBM", SGBM_KWARGS

    if method == "genie3":
        return "RF", RF_KWARGS

    raise ValueError(f"Unsupported method: {method}")


def get_output_separator(output_file):
    """Determine the output separator from the file suffix."""
    if output_file == "-":
        return "\t"

    suffixes = PurePath(output_file).suffixes
    return suffixes_to_separator(suffixes)


def main():
    parser = create_argument_parser()
    args = parser.parse_args()

    if args.num_workers < 1:
        parser.error("--num_workers must be greater than zero.")

    method_params = get_method_parameters(args.method)

    start_time = time.time()

    expression_data = load_exp_matrix(
        args.expression_mtx_fname,
        args.transpose,
        args.sparse,
        args.cell_id_attribute,
        args.gene_attribute,
    )

    if args.sparse:
        gene_names = expression_data[1]
        expression_matrix = expression_data[0]
    else:
        expression_matrix = expression_data
        gene_names = expression_matrix.columns

    print(
        f"Loaded expression matrix with {expression_matrix.shape[0]} cells "
        f"and {expression_matrix.shape[1]} genes in "
        f"{time.time() - start_time:.2f} seconds.",
        file=sys.stderr,
    )

    tf_names = load_tf_names(args.tfs_fname)

    print(f"Loaded {len(tf_names)} transcription factors.", file=sys.stderr)

    expression_matrix, gene_names, tf_names = _prepare_input(
        expression_matrix,
        gene_names,
        tf_names,
    )

    tf_matrix, tf_matrix_gene_names = to_tf_matrix(
        expression_matrix,
        gene_names,
        tf_names,
    )

    print(
        f"Starting {args.method} with {args.num_workers} worker processes.",
        file=sys.stderr,
    )

    start_time = time.time()

    worker_function = partial(
        run_infer_partial_network,
        gene_names=gene_names,
        expression_matrix=expression_matrix,
        tf_matrix=tf_matrix,
        tf_matrix_gene_names=tf_matrix_gene_names,
        method_params=method_params,
        seed=args.seed,
    )

    with Pool(processes=args.num_workers) as pool:
        adjacency_lists = list(
            tqdm(
                pool.imap(
                    worker_function,
                    target_gene_indices(gene_names, target_genes="all"),
                    chunksize=1,
                ),
                total=len(gene_names),
                desc="Inferring regulatory network",
            )
        )

    adjacency_matrix = pd.concat(adjacency_lists).sort_values(
        by="importance",
        ascending=False,
    )

    print(
        f"Network reconstruction completed in "
        f"{time.time() - start_time:.2f} seconds.",
        file=sys.stderr,
    )

    separator = get_output_separator(args.output)

    if args.output == "-":
        adjacency_matrix.to_csv(sys.stdout, index=False, sep=separator)
    else:
        adjacency_matrix.to_csv(args.output, index=False, sep=separator)


if __name__ == "__main__":
    main()
    