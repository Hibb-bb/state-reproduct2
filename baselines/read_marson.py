import anndata as ad
import numpy as np
import scanpy 


adata = ad.read_h5ad("./marson/train/train.h5ad")
print(list(adata.obsm_keys()))

scanpy.pp.normalize_total(adata)
scanpy.pp.log1p(adata)

adata.obsm["X_hvg"] = adata.X 
# adata.write_h5ad("./marson/train/train.h5ad")

scanpy.write("./marson/train/train.h5ad", adata)


# adata = ad.read_h5ad("./marson/test/test.h5ad")
# scanpy.pp.normalize_total(adata)
# scanpy.pp.log1p(adata)

# scanpy.write("./marson/test/test.h5ad", adata)


# X = adata.X
# if hasattr(X, "toarray"):
#     X = X.toarray()
# X = np.asarray(X)

# adata["obsm"]["X"] = X # put things to obsm
# scanpy.write("./marson/train/train.h5ad", adata)


# print("X.dtype:", X.dtype)
# print("X.shape:", X.shape)
# print("X.min():", X.min())
# print("X.max():", X.max())
# print("X.mean():", X.mean())
# print("Sample of first cell (first 10 genes):", X[0, :10])