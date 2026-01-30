from huggingface_hub import snapshot_download
import os

repo_id = "arcinstitute/Replogle-Nadig-Preprint"
local_dir = "./replogle"
os.makedirs(local_dir, exist_ok=True)

# Download ONLY the specified file
snapshot_download(
    repo_id=repo_id,
    repo_type="dataset",
    local_dir=local_dir,
    local_dir_use_symlinks=False,          # optional: copy files instead of symlinks
    allow_patterns=["replogle.h5ad"],  # restrict to just this file
)