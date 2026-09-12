"""Courteous GPU reservation.

The lab shares these devices. Policy: never take the busiest device, never take more
memory than we declare, and take nothing at all when the stage does not need it.

The oracle stage is pure CPU (exact integer dynamic programming), so it calls nothing
here. This module exists so that later stages reserve deliberately rather than by
default.
"""
from __future__ import annotations

import os
import subprocess


def device_status():
    """(index, total_MiB, used_MiB, util_pct) per visible device, via nvidia-smi."""
    try:
        out = subprocess.check_output(
            ["nvidia-smi",
             "--query-gpu=index,memory.total,memory.used,utilization.gpu",
             "--format=csv,noheader,nounits"],
            stderr=subprocess.DEVNULL).decode()
    except Exception:
        return []
    rows = []
    for line in out.strip().splitlines():
        parts = [p.strip() for p in line.split(",")]
        if len(parts) == 4:
            rows.append(tuple(int(p) for p in parts))
    return rows


def pick_least_used(min_free_mib=2048):
    """Index of the device with the most free memory, or None if none qualifies."""
    rows = device_status()
    if not rows:
        return None
    free = [(total - used, idx) for idx, total, used, _ in rows]
    free.sort(reverse=True)
    best_free, best_idx = free[0]
    return best_idx if best_free >= min_free_mib else None


def reserve(fraction=0.15, min_free_mib=2048, verbose=True):
    """Pin this process to one device and cap its share of that device.

    fraction is a hard cap on the process's allocation, not a request. Returns the
    device index taken, or None if we declined to take one.
    """
    idx = pick_least_used(min_free_mib)
    if idx is None:
        if verbose:
            print("[gpu] no device with enough free memory; staying on CPU")
        return None
    os.environ["CUDA_VISIBLE_DEVICES"] = str(idx)
    try:
        import torch
        if torch.cuda.is_available():
            torch.cuda.set_per_process_memory_fraction(fraction, 0)
            if verbose:
                total = torch.cuda.get_device_properties(0).total_memory / 2**20
                print(f"[gpu] pinned to physical device {idx}, "
                      f"capped at {fraction:.0%} of {total:.0f} MiB")
    except Exception as exc:
        if verbose:
            print(f"[gpu] reservation not applied ({type(exc).__name__}); CPU only")
        return None
    return idx
