#!/usr/bin/env python
"""Simple UMAP visualization to compare adata_real and adata_pred."""

import sys
import argparse
import numpy as np
import scanpy as sc
import anndata as ad
from pathlib import Path

def main():
    parser = argparse.ArgumentParser(description="Plot UMAP comparing real vs predicted data")
    parser.add_argument("--output_dir", type=str, required=True, help="Directory containing adata_real.h5ad and adata_pred.h5ad")
    parser.add_argument("--out_path", type=str, default=None, help="Output path for UMAP plot (default: output_dir/umap_comparison.png)")
    args = parser.parse_args()
    
    output_dir = Path(args.output_dir)
    adata_real_path = output_dir / "adata_real.h5ad"
    adata_pred_path = output_dir / "adata_pred.h5ad"
    
    if not adata_real_path.exists():
        raise FileNotFoundError(f"Could not find {adata_real_path}")
    if not adata_pred_path.exists():
        raise FileNotFoundError(f"Could not find {adata_pred_path}")
    
    # Load data
    print(f"Loading {adata_real_path}...")
    adata_real = sc.read_h5ad(adata_real_path)
    print(f"Loading {adata_pred_path}...")
    adata_pred = sc.read_h5ad(adata_pred_path)
    
    # Sample 2500 cells from each to keep it balanced
    n_sample = 2500
    if len(adata_real) > n_sample:
        print(f"Sampling {n_sample} cells from real (total: {len(adata_real)})...")
        np.random.seed(42)
        idx_real = np.random.choice(len(adata_real), n_sample, replace=False)
        adata_real = adata_real[idx_real].copy()
    if len(adata_pred) > n_sample:
        print(f"Sampling {n_sample} cells from pred (total: {len(adata_pred)})...")
        np.random.seed(42)
        idx_pred = np.random.choice(len(adata_pred), n_sample, replace=False)
        adata_pred = adata_pred[idx_pred].copy()
    
    # Combine with source labels
    adata_real.obs["source"] = "real"
    adata_pred.obs["source"] = "pred"
    combined = ad.concat({"real": adata_real, "pred": adata_pred}, label="source", index_unique="-")
    
    # Compute UMAP
    print("Computing PCA...")
    sc.pp.pca(combined, n_comps=50)
    print("Computing neighbors...")
    sc.pp.neighbors(combined, n_neighbors=15)
    print("Computing UMAP...")
    sc.tl.umap(combined)
    
    # Plot
    out_path = args.out_path or (output_dir / "umap_comparison.png")
    print(f"Saving UMAP plot to {out_path}...")
    fig = sc.pl.umap(combined, color="source", title="UMAP: Real vs Predicted", show=False, return_fig=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    print(f"Done! Saved to {out_path}")

if __name__ == "__main__":
    main()

