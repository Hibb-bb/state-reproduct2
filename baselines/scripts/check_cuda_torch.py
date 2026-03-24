#!/usr/bin/env python3
"""
Fail fast: verify PyTorch can initialize CUDA against the installed NVIDIA driver.
Run before heavy data preprocessing. Set SKIP_CUDA_CHECK=1 to skip (e.g. CPU-only debug).
"""
import os
import sys


def main() -> int:
    if os.environ.get("SKIP_CUDA_CHECK", "").strip() in ("1", "true", "yes"):
        print("[check_cuda_torch] SKIP_CUDA_CHECK set — skipping GPU check.")
        return 0

    try:
        import torch
    except ImportError as e:
        print(f"[check_cuda_torch] ERROR: cannot import torch: {e}", file=sys.stderr)
        return 1

    ver = torch.__version__
    cuda_build = getattr(torch.version, "cuda", None) or "None"
    print(f"[check_cuda_torch] torch {ver}, cuda (build) {cuda_build}")

    if not torch.cuda.is_available():
        print(
            "[check_cuda_torch] ERROR: torch.cuda.is_available() is False. "
            "Install a CUDA build of PyTorch or fix your driver.",
            file=sys.stderr,
        )
        return 1

    try:
        torch.cuda.init()
        _ = torch.cuda.current_device()
        name = torch.cuda.get_device_name(0)
        major, minor = torch.cuda.get_device_capability(0)
    except RuntimeError as e:
        print(
            "[check_cuda_torch] ERROR: driver / PyTorch CUDA mismatch (typical fix: "
            "matching torch wheel for your driver, or upgrade NVIDIA driver):\n"
            f"  {e}",
            file=sys.stderr,
        )
        return 1

    print(f"[check_cuda_torch] OK: GPU 0 {name!r}, capability sm_{major}{minor}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
