"""Main functions for the OmniBenchmark module."""

from pathlib import Path
import gzip
import scanpy as sc
import numpy as np
import h5py
from scipy import sparse

def _decode(x):
    return np.array([
        v.decode("utf-8") if isinstance(v, bytes) else str(v)
        for v in x
    ])


def read_tenx_matrix(h5_path):
    """Load TENx HDF5 (genes x cells) as AnnData.

    H5 stores genes x cells.
    AnnData stores cells x genes.
    """
    with h5py.File(h5_path, "r") as h5:
        g = h5["matrix"]

        data = g["data"][:]
        indices = g["indices"][:]
        indptr = g["indptr"][:]
        shape = tuple(g["shape"][:])

        if "features" in g and "id" in g["features"]:
            gene_ids = g["features/id"][:]
        elif "genes" in g:
            gene_ids = g["genes"][:]
        else:
            gene_ids = np.array([f"gene_{i}".encode() for i in range(shape[0])])

        if "barcodes" in g:
            cell_ids = g["barcodes"][:]
        else:
            cell_ids = np.array([f"cell_{i}".encode() for i in range(shape[1])])

    gene_ids = _decode(gene_ids)
    cell_ids = _decode(cell_ids)

    m = sparse.csc_matrix((data, indices, indptr), shape=shape)

    adata = sc.AnnData(X=m.T.tocsr())
    adata.obs_names = cell_ids
    adata.var_names = gene_ids
    return adata

def write_tenx_matrix(adata, h5_path):
    """Write AnnData cells x genes as TENx-like HDF5 matrix genes x cells."""
    X = adata.X
    if sparse.issparse(X):
        m = X.T.tocsc()
    else:
        m = sparse.csc_matrix(X.T)

    gene_ids = np.asarray(adata.var_names.astype(str))
    cell_ids = np.asarray(adata.obs_names.astype(str))

    str_dtype = h5py.string_dtype(encoding="utf-8")

    with h5py.File(h5_path, "w") as h5:
        g = h5.create_group("matrix")
        g.create_dataset("data", data=m.data, compression="gzip")
        g.create_dataset("indices", data=m.indices, compression="gzip")
        g.create_dataset("indptr", data=m.indptr, compression="gzip")
        g.create_dataset("shape", data=np.asarray(m.shape, dtype=np.int64))
        g.create_dataset("genes", data=gene_ids.astype(str_dtype))
        g.create_dataset("barcodes", data=cell_ids.astype(str_dtype))

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
def select_by_scanpy_pearson_residuals(adata, number_selected):
    """Select HVGs by Scanpy analytic Pearson residuals."""
    # https://scanpy.readthedocs.io/en/stable/tutorials/experimental/pearson_residuals.html
    adata = adata.copy()

    sc.experimental.pp.highly_variable_genes(adata, flavor="pearson_residuals", n_top_genes=number_selected, clip=None)

    return adata.var_names[adata.var["highly_variable"]].tolist()


def process_data(args):
    """Process data using parsed command-line arguments.

    Args:
        args: Parsed arguments from argparse containing:
            - output_dir: Output directory path
            - name: Module name
            - normalized_h5: Input files for normalized.h5 (CLI: --normalized.h5)
            - rawdata_h5ad: Input files for rawdata.h5ad (CLI: --rawdata.h5ad)

    Note: Input IDs with dots (e.g., 'data.raw') are converted to underscores
          in Python variable names (e.g., 'data_raw') but preserve dots in CLI args.
    """
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Processing module: {args.name}")

    # Access stage inputs
    normalized_h5_files = args.normalized_h5[0]
    rawdata_h5ad_files = args.rawdata_h5ad[0]
    number_selected = int(args.number_selected)
    print(f"  rawdata.h5ad: {rawdata_h5ad_files}")
    print(f"  normalized.h5: {normalized_h5_files}")
    print(f"  selection_type: {args.selection_type}")
    print(f"  number_selected: {number_selected}")

    adata_norm = read_tenx_matrix(normalized_h5_files)

    if number_selected > adata_norm.n_vars:
        raise ValueError(
            f"number_selected={number_selected} is larger than number of features={adata_norm.n_vars}"
        )
    
    if args.selection_type == "scanpy_hvg":
        sel_feats = select_by_scanpy_hvg(adata_norm, number_selected)

    # TODO：order by gini coef and select top  N; currently it is based on pvalue
    # elif args.selection_type == "giniclust3":
    #     sel_feats = select_by_giniclust3(adata, number_selected)

    elif args.selection_type == "pearson_residuals":
        adata_raw = sc.read_h5ad(rawdata_h5ad_files)
        shared_cells = adata_norm.obs_names.intersection(adata_raw.obs_names)
        shared_genes = adata_norm.var_names.intersection(adata_raw.var_names)
        adata_filtered = adata_raw[adata_norm.obs_names, adata_norm.var_names].copy()
        sel_feats = select_by_scanpy_pearson_residuals(adata_filtered, number_selected)

    else:
        raise ValueError(f"Unknown selection_type: {args.selection_type}")

    print(f"length(sel_feats): {len(sel_feats)}")

    # Write a simple output file
    adata_selected = adata_norm[:, sel_feats].copy()

    output_file = output_dir / f"{args.name}_normalized_selected.h5"
    print(f"output_file: {output_file}")

    write_tenx_matrix(adata_selected, output_file)

    stat = output_file.stat()
    print(f"Results written to: {output_file}")
    print(f"size: {stat.st_size}")
