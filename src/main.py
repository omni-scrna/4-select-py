"""Main functions for the OmniBenchmark module."""

from pathlib import Path
import gzip
import scanpy as sc
# import giniclust3 as gc

# def select_by_giniclust3(adata, number_selected):
#     """Select features with highest Gini coefficient."""
#     # https://github.com/rdong08/GiniClust3
#     adata = adata.copy()
#     # GiniClust3 uses p-value and min_gini_value cutoffs.
#     # The output size may differ from number_selected.
#     gc.gini.calGini(adata)

#     if "gini" not in adata.var:
#         raise ValueError(
#             f"GiniClust3 did not produce 'gini'. Columns: {adata.var.columns.tolist()}"
#         )

#     return adata.var_names[adata.var["gini"]].tolist()


def select_by_scanpy_hvg(adata, number_selected):
    """Select HVGs using Scanpy's Seurat-like normalized-data method."""
    adata = adata.copy()
    sc.pp.highly_variable_genes(adata, n_top_genes=number_selected, flavor="seurat")

    selected = adata.var_names[adata.var["highly_variable"]].tolist()
    return selected[:number_selected]

# This method expects raw counts as input.
# def select_by_scanpy_pearson_residuals(adata, number_selected):
#     """Select HVGs by Scanpy analytic Pearson residuals."""
#     # https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html
#     adata = adata.copy()

#     sc.experimental.pp.highly_variable_genes(adata, flavor="pearson_residuals", n_top_genes=number_selected, clip=None)

#     return adata.var_names[adata.var["highly_variable"]].tolist()


def process_data(args):
    """Process data using parsed command-line arguments.

    Args:
        args: Parsed arguments from argparse containing:
            - output_dir: Output directory path
            - name: Module name
            - normalized_h5: Input files for normalized.h5 (CLI: --normalized.h5)

    Note: Input IDs with dots (e.g., 'data.raw') are converted to underscores
          in Python variable names (e.g., 'data_raw') but preserve dots in CLI args.
    """
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing module: {args.name}")

    # Access stage inputs
    normalized_h5_files = args.normalized_h5[0]
    number_selected = int(args.number_selected)
    print(f"  normalized.h5: {normalized_h5_files}")
    print(f"  selection_type: {args.selection_type}")
    print(f"  number_selected: {number_selected}")

    adata = sc.read_10x_h5(normalized_h5_files, gex_only=False)
    # TODO: Implement your processing logic here
    # Example: Read inputs, process, write outputs
    if number_selected > adata.n_vars:
        raise ValueError(
            f"number_selected={number_selected} is larger than number of features={adata.n_vars}"
        )
    
    if args.selection_type == "scanpy_hvg":
        sel_feats = select_by_scanpy_hvg(adata, number_selected)

    # TODO：order by gini coef and select top  N; currently it is based on pvalue
    # elif args.selection_type == "giniclust3":
    #     sel_feats = select_by_giniclust3(adata, number_selected)

    # TODO: input supposed to be the raw counts
    # elif args.selection_type == "pearson_residuals":
    #     sel_feats = select_by_scanpy_pearson_residuals(adata, number_selected)

    else:
        raise ValueError(f"Unknown selection_type: {args.selection_type}")

    print(f"length(sel_feats): {len(sel_feats)}")

    # Write a simple output file
    output_file = output_dir / f"{args.name}_selected.txt.gz"
    with gzip.open(output_file, "wt") as f:
        for feat in sel_feats:
            f.write(f"{feat}\n")

    print(f"Results written to: {output_file}")
