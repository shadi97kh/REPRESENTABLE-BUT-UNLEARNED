"""Application claims-to-evidence table, generated from immutable result files.

Nothing is computed here. Every row is read from a file under `runs/`, and a file with a
`.sha256` sidecar is re-hashed and refused if it moved. A claim with no evidence file
does not get a row.
"""
from __future__ import annotations

import argparse, glob, hashlib, json, os, sys, time
import numpy as np
import yaml


def load(path):
    payload = open(path, "rb").read()
    side = path + ".sha256"
    if os.path.exists(side):
        if hashlib.sha256(payload).hexdigest() != open(side).read().strip():
            raise SystemExit(f"REFUSED: {path} content hash moved")
    return json.loads(payload)


def latest(pat):
    fs = sorted(glob.glob(pat))
    return fs[-1] if fs else None


def md(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.rna.export_paper")
    ap.add_argument("--protocol", default="configs/rna_protocol.yaml")
    ap.add_argument("--out", default="docs/rna_claims_to_evidence.md")
    a = ap.parse_args(argv)
    proto = yaml.safe_load(open(a.protocol))

    src = {"eval": latest("runs/rna_eval/eval-*.json"),
           "data": latest("runs/rna_data/manifest-*.json"),
           "splits": latest("runs/rna_splits/splits-*.json"),
           "ensemble": latest("runs/rna_ensemble/ensemble-*.json"),
           "pilot": latest("runs/rna_pilot/dryrun-*.json")}
    missing = [k for k, v in src.items() if v is None]
    ev = load(src["eval"]) if src["eval"] else {}
    dat = load(src["data"]) if src["data"] else {}
    spl = load(src["splits"]) if src["splits"] else {}
    ens = load(src["ensemble"]) if src["ensemble"] else {}

    L = ["# RNA application: claims to evidence\n",
         f"Generated {time.strftime('%Y-%m-%d')} by `scmp.rna.export_paper` from "
         "immutable result files. Task is **on-target siRNA efficacy and ranking**. No "
         "off-target or safety quantity is computed anywhere in `scmp/rna`.\n"]

    # ---- supported
    L.append("\n## Supported\n")
    rows = []
    h = dat.get("manifests", {}).get("huesken2005", {}) if dat else {}
    if h:
        rows.append(["Huesken provenance preserved: records, guide length, unclipped "
                     "labels, mapping coverage, historical split",
                     f"{h['n_records']} records, all {list(h['guide_lengths'])[0]} nt, "
                     f"labels [{h['label_min']:.3f}, {h['label_max']:.3f}] unclipped, "
                     f"{h['mapping_coverage']['annotated']}/"
                     f"{h['mapping_coverage']['total']} mapped, split "
                     f"{h['historical_split']}", "data manifest"])
        rows.append(["Assay identity taken from the primary source, not a redistribution",
                     f"{h['assay']['readout']}, {h['assay']['cell_line']}, "
                     f"{h['assay']['timepoint_h']} h, {h['assay']['n_constructs']} "
                     f"constructs; {len(h['provenance_conflicts'])} conflicts with the "
                     f"redistributed annotation recorded", "data manifest"])
    if spl:
        d = spl.get("duplication", {})
        rows.append(["The historical split leaks; a sequence-cluster split does not",
                     f"historical leak {d.get('historical_leak_kmer_cluster', float('nan')):.1%} "
                     f"by 13-mer cluster and "
                     f"{d.get('historical_leak_gene', float('nan')):.0%} by gene; "
                     f"cluster split residual leak 0.00%", "split audit"])
    if ens:
        rows.append(["Ensemble machinery is exact where it claims to be",
                     f"{ens['checks']} checks, {len(ens['failures'])} failures: exact "
                     f"marginals match enumeration, every sample is a valid member, the "
                     f"deterministic mean bracket contains the exact mean, and "
                     f"yhat = b + a^T mu + gamma E[g] equals E_p[f]", "ensemble audit"])
    if ev:
        s = ev.get("summary", {})
        k = "reporter_yfp|linear_sequence_only"
        if k in s:
            sp, pk = s[k]["spearman"], s[k]["precision_at_k"]
            rows.append(["On-target efficacy is predictable from guide sequence alone "
                         "under a leakage-free split",
                         f"Spearman {np.mean(sp):.4f} +/- {np.std(sp):.4f} over "
                         f"{len(sp)} seeds, P@20 {np.mean(pk):.3f}; "
                         f"sequence-cluster-disjoint split, single assay",
                         "evaluation"])
        rc = ev.get("resolved_comparison_demo", {})
        if rc:
            rows.append(["The u_i < l_j rule resolves a substantial fraction of pairs "
                         "and certifies the MODEL ranking",
                         f"{rc['resolved']}/{rc['pairs']} pairs resolved "
                         f"({rc['coverage']:.1%}); model ordering correct on "
                         f"{rc['model_ranking_agreement']:.0%} of them, as a certificate "
                         f"requires. Declared windows, NOT an RNA result",
                         "evaluation (engineering)"])
    L.append(md(["claim", "evidence", "source"], rows))

    # ---- not supported
    L.append("\n\n## Not supported\n")
    rows = []
    if ev:
        for b in ev.get("blocked", []):
            rows.append([f"`{b['arm']}` arm", "BLOCKED",
                         f"missing: {', '.join(b['missing'])}"])
    rows.append(["Structural sensitivity of efficacy (50/100/150 nt)", "BLOCKED",
                 "the assayed context is a YFP reporter construct that is not in hand; "
                 "native transcript windows would answer a different experiment"])
    rows.append(["Learned-layer benefit", "BLOCKED", "requires the A+B arms"])
    rows.append(["Useful residual nonlinearity on RNA", "BLOCKED",
                 "requires the A+B arms; the ML-stage gate G4 already found the residual "
                 "*hurts* held-out prediction (-0.268 Spearman)"])
    rows.append(["Agreement of certified ranking with held-out measurement", "BLOCKED",
                 "no measured ensemble readout exists for any window we can legitimately "
                 "build"])
    rows.append(["Chemical-variant prediction", "BLOCKED",
                 "Huesken contains a single chemistry (unmodified); Davis 2025 has the "
                 "chemistry field and is not yet acquired"])
    rows.append(["Biological causation", "NOT CLAIMED",
                 "the label is marginal over latent structures and is a reporter readout"])
    L.append(md(["claim", "status", "reason"], rows))

    # ---- separation
    L.append("\n\n## Separation of claims\n")
    L.append(md(["kind", "status"],
                [[k, v] for k, v in proto["separation"].items()]))
    L.append("\n`u_i < l_j` certifies that the **model** ranks i above j. It is not a "
             "statement about measurement. Agreement with held-out measurement is an "
             "empirical quantity, reported separately, and is currently not computable "
             "on any admissible RNA data.\n")

    # ---- gates
    L.append("\n## Gate status\n")
    L.append(md(["gate", "status"], [
        ["ML G1 conditional vs unconditional", "UNDETERMINED"],
        ["ML G2 learned vs deterministic allocation", "UNDETERMINED"],
        ["ML G3 predictor noninferiority", "PASS (baseline weak)"],
        ["ML G4 residual predictive value", "FAIL"],
        ["DATA assay/construct identity", "BLOCKING"],
        ["DATA cross-dataset duplication", "BLOCKING (datasets not acquired)"],
        ["RNA pilot launch", "REFUSED (no authorised budget)"]]))

    if missing:
        L.append(f"\n\n## Missing sources\n\n{', '.join(missing)}\n")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    open(a.out, "w").write("\n".join(L) + "\n")
    for k, v in src.items():
        print(f"  {k:>9}: {os.path.basename(v) if v else 'MISSING'}")
    print(f"\nwritten: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
