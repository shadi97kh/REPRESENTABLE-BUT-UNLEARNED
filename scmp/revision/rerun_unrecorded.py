"""Create immutable records for audits that persist nothing on success.

`audit_bounds`, `audit_conditional` and `check_certificates` write a JSON record
only when they FIND A VIOLATION. A clean run prints to stdout and leaves no
artefact, so the run report's headline soundness numbers -- 1,262,952 bound
checks, 72/72 certificates -- cannot be traced to any file.

This module does NOT modify those audits. It runs them unchanged, captures
stdout, exit status, timing and the source hashes of the modules involved, and
writes that as a new record under `runs/revision/`. The existing code stays
exactly as it produced the historical results; the missing artefact is added
alongside rather than by editing history.

Parsed fields are marked `parsed_from_stdout` so nobody mistakes a scraped
number for a structured result.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
import time

from .provenance import sha256_file

JOBS = {
    "audit_bounds": {
        "cmd": [sys.executable, "-m", "scmp.audit_bounds",
                "--config", "configs/exhaustive.yaml"],
        "sources": ["scmp/audit_bounds.py", "scmp/bounds.py", "scmp/model.py"],
        "config": "configs/exhaustive.yaml",
        "patterns": {
            "cases": r"cases\s*:\s*(\d+)",
            "feasible_structures": r"\((\d+) feasible structures visited\)",
            "intermediate_bracket_checks": r"intermediate bracket checks:\s*(\d+)",
            "reference_forward_checks": r"reference-forward checks\s*:\s*(\d+)",
            "masked_bound_checks": r"masked-bound checks\s*:\s*(\d+)",
            "gamma0_exact": r"gamma=0 cases exact\s*:\s*(\d+)",
        },
        "flags": {"no_violation": "no soundness violation found"},
    },
    "check_certificates": {
        "cmd": [sys.executable, "-m", "scmp.check_certificates",
                "--run", "runs/tiny_complete"],
        "sources": ["scmp/check_certificates.py", "scmp/refine.py"],
        "config": "runs/tiny_complete",
        "patterns": {
            "certificates_checked": r"certificate check:\s*(\d+) file",
            "obligations_discharged": r"obligations discharged\s*:\s*(\d+)",
            "accepted": r"certificates accepted\s*:\s*(\d+)/",
            "accepted_of": r"certificates accepted\s*:\s*\d+/(\d+)",
        },
        "flags": {"all_accepted": "all certificates accepted"},
    },
    "audit_conditional": {
        "cmd": [sys.executable, "-m", "scmp.audit_conditional"],
        "sources": ["scmp/audit_conditional.py", "scmp/conditional.py"],
        "config": "configs/exhaustive.yaml",
        "patterns": {
            "cases": r"cases audited\s*:\s*(\d+)",
            "feasible_structures": r"\((\d+) feasible structures\)",
            "oracle_calls": r"oracle calls\s*:\s*(\d+)",
            "tighter": r"tighter than baseline\s*:\s*(\d+)",
            "proved_impossible": r"edges proved impossible\s*:\s*(\d+)",
            "fallback_wins": r"fallback won at top level\s*:\s*(\d+)",
        },
        "flags": {"no_violation": "no soundness violation found"},
    },
}


def run_job(name, spec, timeout):
    t0 = time.time()
    proc = subprocess.run(spec["cmd"], capture_output=True, text=True,
                          timeout=timeout)
    out = proc.stdout
    parsed = {}
    for key, pat in spec["patterns"].items():
        m = re.search(pat, out, re.IGNORECASE)
        parsed[key] = int(m.group(1).replace(",", "")) if m else None
    flags = {k: (v in out) for k, v in spec.get("flags", {}).items()}
    unparsed = [k for k, v in parsed.items() if v is None]
    return {"job": name, "command": " ".join(spec["cmd"]), "flags": flags,
            "unparsed_keys": unparsed,
            "returncode": proc.returncode, "seconds": round(time.time() - t0, 2),
            "config": spec["config"],
            "source_hashes": {p: sha256_file(p) for p in spec["sources"]
                              if os.path.exists(p)},
            "parsed_from_stdout": parsed,
            "stdout": out, "stderr_tail": proc.stderr[-2000:],
            "note": "the audit was run UNCHANGED; this record is added alongside "
                    "it. Parsed values are scraped from stdout and are weaker "
                    "evidence than a structured result would be."}


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.rerun_unrecorded",
        description="Run clean-on-success audits and persist their output as an "
                    "immutable record. Does not modify the audits.")
    ap.add_argument("--jobs", nargs="*", default=sorted(JOBS),
                    choices=sorted(JOBS))
    ap.add_argument("--timeout", type=int, default=3600)
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--reparse", metavar="RECORD",
                    help="re-parse a previous record's saved stdout instead of "
                         "re-executing the audits")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)

    if a.dry_run:
        print("dry run -- would execute, CPU only, no GPU:")
        for j in a.jobs:
            print(f"  {j:<20} {' '.join(JOBS[j]['cmd'])}")
        print(f"\nbudget: CPU-only, single process each, timeout {a.timeout}s per "
              f"job. No GPU hours consumed. No historical record is modified.")
        return 0

    if a.reparse:
        prev = json.load(open(a.reparse))
        out_rows = []
        for r in prev["results"]:
            spec = JOBS[r["job"]]
            parsed = {}
            for key, pat in spec["patterns"].items():
                m = re.search(pat, r.get("stdout", ""), re.IGNORECASE)
                parsed[key] = int(m.group(1).replace(",", "")) if m else None
            r["parsed_from_stdout"] = parsed
            r["flags"] = {k: (v in r.get("stdout", ""))
                          for k, v in spec.get("flags", {}).items()}
            r["unparsed_keys"] = [k for k, v in parsed.items() if v is None]
            r["reparsed_from"] = a.reparse
            out_rows.append(r)
            print(f"{r['job']}: {parsed}  flags={r['flags']}")
        os.makedirs(a.out_dir, exist_ok=True)
        path = os.path.join(
            a.out_dir, f"rerun_unrecorded-{time.strftime('%Y%m%d-%H%M%S')}.json")
        with open(path, "w") as fh:
            json.dump({"kind": "revision_rerun_unrecorded",
                       "created": time.strftime("%Y-%m-%d %H:%M:%S"),
                       "reparsed_from": a.reparse, "gpu_used": False,
                       "cpu_seconds_total": 0.0, "results": out_rows}, fh, indent=1)
        with open(path + ".sha256", "w") as fh:
            fh.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
        print(f"\nwritten {path}")
        return 0

    results, t0 = [], time.time()
    for j in a.jobs:
        print(f"running {j} ...", flush=True)
        try:
            r = run_job(j, JOBS[j], a.timeout)
        except subprocess.TimeoutExpired:
            r = {"job": j, "returncode": None, "timeout_s": a.timeout,
                 "note": "timed out; no record of a clean run created"}
        results.append(r)
        print(f"  exit {r.get('returncode')}  {r.get('seconds')}s  "
              f"{r.get('parsed_from_stdout')}")

    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir,
                        f"rerun_unrecorded-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump({"kind": "revision_rerun_unrecorded",
                   "created": time.strftime("%Y-%m-%d %H:%M:%S"),
                   "cpu_seconds_total": round(time.time() - t0, 2),
                   "gpu_used": False, "results": results}, fh, indent=1)
    with open(path + ".sha256", "w") as fh:
        fh.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
    print(f"\nwritten {path}")
    return 0 if all(r.get("returncode") == 0 for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
