"""Snapshot the exact state a revision run was made from.

The phase-one records could not distinguish two F7 runs that shared a git SHA and
a config hash but differed materially, because the tree was dirty and the git SHA
says nothing about uncommitted edits. This module fixes that failure mode: it
hashes the CONTENT of every source, config and data file actually present,
whether or not it is committed, and records the dirty flag explicitly.

A snapshot is written under `runs/revision/snapshot-*.json` and is referenced by
id from every revision result, so a later reader can tell exactly which bytes
produced a number even if the tree has since moved on.
"""
from __future__ import annotations

import hashlib
import importlib.metadata as md
import json
import os
import platform
import subprocess
import sys
import time

SOURCE_DIRS = ("certmp", "scmp", "experiments", "tests")
CONFIG_DIRS = ("configs",)
DATA_DIRS = ("data",)
SKIP_DIRS = {"__pycache__", ".git", ".pytest_cache", ".ipynb_checkpoints"}
DATA_MAX_BYTES = 64 << 20      # hash data files up to 64 MB; record size beyond


def sha256_file(path: str, n: int = 64) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()[:n]


def _walk(root: str, suffixes: tuple | None = None):
    if not os.path.isdir(root):
        return
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = sorted(d for d in dirnames if d not in SKIP_DIRS)
        for fn in sorted(filenames):
            if suffixes and not fn.endswith(suffixes):
                continue
            yield os.path.join(dirpath, fn)


def _git(*args) -> str:
    try:
        return subprocess.run(["git", *args], capture_output=True, text=True,
                              timeout=30).stdout.strip()
    except Exception:
        return ""


def git_state() -> dict:
    """Commit, dirty flag, and the per-path status of everything not committed.

    `dirty` being true is not a defect; failing to record WHICH files were dirty
    is. The porcelain listing is kept verbatim so an uncommitted edit is
    attributable to a file rather than to the tree as a whole.
    """
    porcelain = _git("status", "--porcelain")
    entries = [{"status": ln[:2].strip(), "path": ln[3:]}
               for ln in porcelain.splitlines() if ln.strip()]
    return {"commit": _git("rev-parse", "HEAD"),
            "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
            "dirty": bool(entries),
            "uncommitted": entries,
            "note": "a commit SHA alone cannot identify a dirty tree; the file "
                    "hashes below are the authoritative record"}


def environment() -> dict:
    pkgs = {}
    for name in ("numpy", "torch", "scipy", "ViennaRNA", "pyyaml", "reportlab"):
        try:
            pkgs[name] = md.version(name)
        except Exception:
            pkgs[name] = None
    try:
        import torch
        cuda = {"available": bool(torch.cuda.is_available()),
                "device_count": int(torch.cuda.device_count())
                if torch.cuda.is_available() else 0,
                "devices": [torch.cuda.get_device_name(i)
                            for i in range(torch.cuda.device_count())]
                if torch.cuda.is_available() else []}
    except Exception:
        cuda = {"available": False, "device_count": 0, "devices": []}
    return {"python": sys.version.split()[0], "platform": platform.platform(),
            "machine": platform.machine(), "cpu_count": os.cpu_count(),
            "packages": pkgs, "cuda": cuda,
            "env_vars": {k: os.environ.get(k) for k in
                         ("OMP_NUM_THREADS", "MKL_NUM_THREADS",
                          "CUDA_VISIBLE_DEVICES", "PYTHONHASHSEED")}}


def file_hashes() -> dict:
    """Content hashes for source, config and data. Uncommitted files included."""
    out = {"source": [], "config": [], "data": []}
    for d in SOURCE_DIRS:
        for p in _walk(d, (".py",)):
            out["source"].append({"path": p, "sha256": sha256_file(p),
                                  "bytes": os.path.getsize(p)})
    for d in CONFIG_DIRS:
        for p in _walk(d, (".yaml", ".yml", ".json")):
            out["config"].append({"path": p, "sha256": sha256_file(p),
                                  "bytes": os.path.getsize(p)})
    for d in DATA_DIRS:
        for p in _walk(d):
            size = os.path.getsize(p)
            rec = {"path": p, "bytes": size}
            if size <= DATA_MAX_BYTES:
                rec["sha256"] = sha256_file(p)
            else:
                rec["sha256"] = None
                rec["note"] = f"larger than {DATA_MAX_BYTES} bytes; not hashed"
            out["data"].append(rec)
    return out


def snapshot(seeds=None, note: str = "", out_dir: str = "runs/revision") -> dict:
    """Record the full state and write it. Returns the snapshot including its id."""
    fh = file_hashes()
    # One digest over every content hash, so two snapshots can be compared in O(1)
    # and a single changed byte anywhere changes the id.
    agg = hashlib.sha256()
    for kind in ("source", "config", "data"):
        for rec in fh[kind]:
            agg.update(f"{kind}:{rec['path']}:{rec['sha256']}".encode())
    snap = {"kind": "revision_snapshot",
            "created": time.strftime("%Y-%m-%d %H:%M:%S"),
            "tree_digest": agg.hexdigest(),
            "git": git_state(),
            "environment": environment(),
            "seeds": list(seeds) if seeds is not None else None,
            "counts": {k: len(v) for k, v in fh.items()},
            "files": fh,
            "note": note}
    os.makedirs(out_dir, exist_ok=True)
    path = os.path.join(out_dir, f"snapshot-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as f:
        json.dump(snap, f, indent=1, sort_keys=True)
    with open(path + ".sha256", "w") as f:
        f.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
    snap["path"] = path
    return snap


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.provenance",
        description="Hash source (including uncommitted), config, data, "
                    "environment and seeds into an immutable snapshot record.")
    ap.add_argument("--seeds", type=int, nargs="*", default=None)
    ap.add_argument("--note", default="")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)
    s = snapshot(a.seeds, a.note, a.out_dir)
    g = s["git"]
    print(f"snapshot      {s['path']}")
    print(f"tree digest   {s['tree_digest'][:32]}")
    print(f"git           {g['commit'][:12]} on {g['branch']}  "
          f"dirty={g['dirty']} ({len(g['uncommitted'])} uncommitted paths)")
    print(f"files         {s['counts']['source']} source, "
          f"{s['counts']['config']} config, {s['counts']['data']} data")
    e = s["environment"]
    print(f"python        {e['python']}  numpy {e['packages'].get('numpy')}  "
          f"torch {e['packages'].get('torch')}")
    print(f"cuda          {e['cuda']['device_count']} device(s) {e['cuda']['devices']}")
    print(f"seeds         {s['seeds']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
