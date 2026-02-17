#!/usr/bin/env python3
"""
Fix Marson h5ad files to only include unique perturbations that actually appear in cells.
This reduces the categories array from 1M+ to ~10k, making it manageable.
"""
import scanpy as sc
import numpy as np
from pathlib import Path

def fix_categories(h5ad_path: Path, output_path: Path, pert_col: str = "guide_target_gene_symbol"):
    """Fix the categories array to only include perturbations that appear in cells."""
    print(f"Processing {h5ad_path}...")
    
    # Read the data
    adata = sc.read_h5ad(h5ad_path)
    
    # Get unique perturbations that actually appear in cells
    unique_perts = adata.obs[pert_col].unique()
    print(f"  Total cells: {len(adata)}")
    print(f"  Unique perturbations in cells: {len(unique_perts)}")
    
    # Check current categories
    if adata.obs[pert_col].dtype.name == 'category':
        old_categories = adata.obs[pert_col].cat.categories
        print(f"  Old categories count: {len(old_categories)}")
        
        # Re-categorize to only include unique values that appear
        adata.obs[pert_col] = adata.obs[pert_col].astype(str)
        adata.obs[pert_col] = adata.obs[pert_col].astype('category')
        
        new_categories = adata.obs[pert_col].cat.categories
        print(f"  New categories count: {len(new_categories)}")
        print(f"  Reduction: {len(old_categories)} -> {len(new_categories)}")
    
    # Create output directory if it doesn't exist
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save the fixed file
    print(f"  Saving to {output_path}...")
    adata.write_h5ad(output_path)
    print(f"  Done!")
    
    return len(unique_perts)

if __name__ == "__main__":
    base_dir = Path("/mnt/data/Marson/ML_splits/temp_data")
    output_dir = base_dir / "fixed"
    
    # Fix train and test files
    train_path = base_dir / "train" / "adata.h5ad"
    test_path = base_dir / "test" / "adata.h5ad"
    
    train_output = output_dir / "train" / "adata.h5ad"
    test_output = output_dir / "test" / "adata.h5ad"
    
    print("Fixing Marson h5ad files to only include unique perturbations in cells...")
    print(f"Output directory: {output_dir}\n")
    
    if train_path.exists():
        fix_categories(train_path, train_output)
        print()
    
    if test_path.exists():
        fix_categories(test_path, test_output)
        print()
    
    print("All files fixed!")
    print(f"\nFixed files saved to: {output_dir}")
    print("Update marson.toml to point to: /mnt/data/Marson/ML_splits/temp_data/fixed")

