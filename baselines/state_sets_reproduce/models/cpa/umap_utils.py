import scanpy as sc
import anndata as ad
from matplotlib import pyplot as plt

def plot_umap(real, pred, out_path, title="UMAP: Real vs Reconstructed"):
    """Simple UMAP plot for normalized data. Real and pred should be (n_cells, n_genes) arrays."""
    combined = ad.concat({"real": ad.AnnData(real), "pred": ad.AnnData(pred)}, label="source", index_unique="-")
    sc.pp.pca(combined, n_comps=50); sc.pp.neighbors(combined, n_neighbors=15); sc.tl.umap(combined)
    fig = sc.pl.umap(combined, color="source", title=title, show=False, return_fig=True)
    fig.savefig(out_path, dpi=150, bbox_inches="tight"); plt.close(fig)

