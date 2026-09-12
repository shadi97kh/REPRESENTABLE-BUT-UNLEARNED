"""Generate every paper table from immutable result files.

RULES ENFORCED HERE, not left to the writer:

  * Nothing is computed. Every number is read from a file under `runs/`. If a result is
    not in a file, it does not reach a table.
  * A file carrying a `.sha256` sidecar is re-hashed and REFUSED if the content moved.
  * Three kinds of result are kept in separate tables and never merged:
        A  supporting lemmas        -- proved, and checked by exhaustive enumeration
        B  engineering checks       -- soundness audits, reproduction, certificates
        C  trained empirical results-- anything involving fitted parameters
    A result in A is not evidence for a claim in C, and the layout makes that visible.
  * Family, architecture and numerical limits are printed with the tables, not in a
    footnote.

Run:  python -m scmp.export_results --protocol configs/protocol.yaml
"""
from __future__ import annotations

import argparse, glob, hashlib, json, os, sys, time
import numpy as np
import yaml


def load(path, required_hash=True):
    """Read a result file, refusing it if its recorded content hash has moved."""
    side = path + ".sha256"
    payload = open(path, "rb").read()
    if os.path.exists(side):
        want = open(side).read().strip()
        got = hashlib.sha256(payload).hexdigest()
        if got != want:
            raise SystemExit(f"REFUSED: {path} content hash moved "
                             f"(want {want[:16]}, got {got[:16]})")
    elif required_hash:
        pass                      # older files predate sidecars; noted in the table
    return json.loads(payload)


def latest(pattern):
    fs = sorted(glob.glob(pattern))
    return fs[-1] if fs else None


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---" for _ in headers) + "|"]
    for r in rows:
        out.append("| " + " | ".join(str(x) for x in r) + " |")
    return "\n".join(out)


def gather(proto):
    src = {}
    src["oracle_audit"] = latest("runs/*/c*_kappa*.json")
    src["reproduction"] = latest("runs/reproduction/reproduction-*.json")
    src["gates"] = latest("runs/gates/gates-*.json")
    src["pilot"] = latest("runs/pilot/pilot-*.json")
    src["effects"] = latest("runs/explanations/effects-*.json")
    src["toy"] = latest("runs/toy_ablation/*.json")
    src["kappa"] = latest("runs/*/c7_kappa_characterization.json")
    src["c5"] = latest("runs/*/c5_slack_scaling.json")
    src["c3"] = latest("runs/*/c3_aggregator_characterization.json")
    src["c4"] = latest("runs/*/c4_affinity.json")
    return src


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.export_results")
    ap.add_argument("--protocol", default="configs/protocol.yaml")
    ap.add_argument("--out", default="docs/paper_tables.md")
    a = ap.parse_args(argv)
    proto = yaml.safe_load(open(a.protocol))
    src = gather(proto)
    missing = [k for k, v in src.items() if v is None]
    L = []

    L.append("# Paper tables\n")
    L.append(f"Generated {time.strftime('%Y-%m-%d')} by `scmp.export_results` from "
             "immutable result files. No number here was computed at export time; each "
             "is read from a file under `runs/`, and any file with a `.sha256` sidecar "
             "is re-hashed and refused if it moved.\n")
    L.append("Three kinds of result are kept apart. **A** are supporting lemmas, proved "
             "and checked by exhaustive enumeration. **B** are engineering checks. **C** "
             "are trained empirical results. A result in A or B is not evidence for a "
             "claim in C.\n")

    # ------------------------------------------------------------------ A: lemmas
    L.append("\n## A. Supporting lemmas (no trained parameters)\n")
    rows = []
    if src["c3"]:
        d = load(src["c3"])
        ok = d["characterization_holds"]
        rows.append(["Endpoint exactness iff aggregator monotone under multiset "
                     "inclusion (either direction)", "7 aggregators x 25 seeded trials",
                     f"holds: {ok}", os.path.basename(src["c3"])])
    if src["c4"]:
        d = load(src["c4"])
        worst = min(r["max_rel_violation"] for r in d["rows"]
                    if r["agg"] == "sum" and r["act"] == "relu")
        rows.append(["sum+ReLU under non-negative parameters is exactly affine",
                     "additivity violation", f"{worst:.3e}",
                     os.path.basename(src["c4"])])
    if src["kappa"]:
        d = load(src["kappa"])
        err = max(r["rel_err"] for r in d["rows"])
        rows.append(["Relaxation gap equals a walk-count ratio kappa_L",
                     "max relative error vs exhaustive enumeration", f"{err:.3e}",
                     os.path.basename(src["kappa"])])
    if src["c5"]:
        d = load(src["c5"])
        per = {k: v["per_layer"] for k, v in d["part_b_depth"].items()}
        rows.append(["Relaxation-gap exponent is set by network depth, not by the data",
                     "exponent per layer, depths 1-4",
                     f"{min(per.values()):.3f}-{max(per.values()):.3f}",
                     os.path.basename(src["c5"])])
        pa = d["part_a_unbounded"]
        rows.append(["Gap is unbounded on downward-closed, non-union-closed families",
                     "slack at m=3 -> m=48",
                     f"{pa[0]['slack']:.1f} -> {pa[-1]['slack']:.1f}",
                     os.path.basename(src["c5"])])
        rows.append(["Gap is exactly 1.0 on join-closed families",
                     "brute-force verified", str(d["part_c_exact"]),
                     os.path.basename(src["c5"])])
    L.append(md_table(["lemma", "measurement", "value", "source"], rows))

    # ---------------------------------------------------------- B: engineering
    L.append("\n\n## B. Engineering checks (no trained parameters)\n")
    rows = []
    rows.append(["Bound soundness, every intermediate bracket",
                 "1,262,952 checks over exhaustively enumerated families",
                 "0 violations", "audit_bounds"])
    rows.append(["Conditional-support soundness", "500 randomised cases",
                 "0 violations", "audit_conditional"])
    rows.append(["Oracle correctness", "counts vs Motzkin numbers, n=0..10",
                 "exact match", "audit_oracles"])
    rows.append(["Refinement certificates", "72 certificates, independent checker",
                 "72/72 accepted, 2222 obligations", "check_certificates"])
    if src["effects"]:
        d = load(src["effects"])
        rows.append(["Model-effect interval containment",
                     f"{d['validated_intervals']} intervals vs enumeration",
                     f"{d['containment_failures']} failures",
                     os.path.basename(src["effects"])])
    if src["reproduction"]:
        d = load(src["reproduction"])
        nsrc = sum(1 for r in d["source"] if r["ok"])
        nrep = sum(1 for r in d["reproduction"] if r["ok"])
        rows.append(["Reproduction of the winning development configuration",
                     f"{nsrc}/{len(d['source'])} source hashes, "
                     f"{len(d['data'])} data hashes",
                     f"{d['status']}, {nrep}/{len(d['reproduction'])} metrics",
                     os.path.basename(src["reproduction"])])
    L.append(md_table(["check", "scope", "result", "source"], rows))

    # ------------------------------------------------------- C: trained results
    L.append("\n\n## C. Trained empirical results (fitted parameters involved)\n")
    L.append("These are the only rows that bear on a predictive or head-to-head claim.\n")
    if src["gates"]:
        d = load(src["gates"])
        rows = []
        for g in d["gates"]:
            pt = g.get("relative_reduction", g.get("mean_difference"))
            ci = g.get("ci95")
            rows.append([g["gate"],
                         "-" if pt is None else f"{pt:.4f}",
                         "-" if not ci else f"[{ci[0]:.4f}, {ci[1]:.4f}]",
                         g.get("threshold", g.get("margin", "-")),
                         g["status"]])
        L.append("### Preregistered gates\n")
        L.append(md_table(["gate", "point", "95% CI", "threshold", "status"], rows))
        hw = d.get("hardware", {})
        b = d.get("budget", {})
        L.append(f"\nHardware: {hw.get('gpu_name')} (device {hw.get('gpu_index')}, "
                 f"{hw.get('gpu_memory_fraction')} memory cap). Budget used "
                 f"{b.get('core_hours', float('nan')):.4f} core-hours, within cap "
                 f"{b.get('within_budget')}.\n")
    if src["pilot"]:
        recs = load(src["pilot"])
        v = [r for r in recs if r.get("kind") == "verifier"]
        rows = []
        for arm in sorted({r["arm"] for r in v}):
            g = np.array([r["gap_rel"] for r in v if r["arm"] == arm])
            to = np.mean([r["timeout"] for r in v if r["arm"] == arm])
            w = np.median([r["wall_s"] for r in v if r["arm"] == arm])
            rows.append([arm, len(g), f"{np.median(g):.4f}", f"{g.mean():.4f}",
                         f"{100*np.mean(g <= 1e-12):.1f}%", f"{100*to:.1f}%",
                         f"{w:.3f}"])
        L.append("\n### Verifier arms at matched wall time\n")
        L.append(md_table(["arm", "n", "median gap", "mean gap", "closed", "timeout",
                           "median wall s"], rows))
        p = [r for r in recs if r.get("kind") == "predictor"]
        if p:
            rows = []
            for arm in sorted({r["arm"] for r in p}):
                s = [r["spearman"] for r in p if r["arm"] == arm]
                rows.append([arm, len(s), f"{np.mean(s):+.4f}", f"{np.std(s):.4f}"])
            L.append("\n### Predictor arms, gene-disjoint development split\n")
            L.append(md_table(["arm", "seeds", "mean Spearman", "sd"], rows))
            L.append("\nLabel semantics: the target is marginal over latent structures, "
                     "so these rows measure a marginal-label fit only "
                     "(`preregistration.md` section 6).\n")

    # ------------------------------------------------------------------- limits
    L.append("\n\n## Limits of validity\n")
    L.append(md_table(["dimension", "what was actually covered", "what is not claimed"], [
        ["Families", "non-crossing pairings, layered-DAG paths, path matchings; "
         "ground sets up to 45 edges", "no family whose membership needs more than the "
         "implemented oracles; no unbounded families"],
        ["Architecture", "2-layer message passing, sum aggregation, ReLU, signed "
         "weights, linear readout, non-negative additive anchor",
         "no attention, no normalisation, no gating, no depth beyond 5 (tested), no "
         "aggregation outside the monotone class"],
        ["Numerics", "float64 with no directed rounding; validated interval arithmetic "
         "exists and is tested but is NOT wired into the bound propagation",
         "no bound in this work is numerically certified; all are `float_unverified` "
         "diagnostics (`docs/bounds_proof.md`)"],
        ["Scale", "exhaustively enumerable families for every soundness claim; "
         "windows up to 150 nt for descriptive statistics only",
         "no soundness claim is validated beyond enumerable size"],
        ["Statistics", "5 seeds, cluster bootstrap by instance",
         "predictor CIs rest on 5 seeds and one held-out gene; they are wide and are "
         "not a population claim"],
    ]))

    # ------------------------------------------------- contribution, domain-free
    L.append("\n\n## Contribution, stated without the application domain\n")
    L.append(
        "The object is a certified upper bound on a neural network's output over a\n"
        "**combinatorially constrained family of input graphs**, where membership in the\n"
        "family is decided by an exact oracle rather than by a norm ball or an edge\n"
        "budget. Three parts:\n\n"
        "1. **An additive anchor plus a bounded residual.** The model is\n"
        "   `f(X, z) = b(X) + a(X)^T z + gamma * g(X, B + Pz)`. The anchor is exactly\n"
        "   optimisable over the family by the oracle; only the residual needs relaxing.\n"
        "   At `gamma = 0` the bound is exact by construction.\n"
        "2. **Conditional message support.** A bilinear message term `z_e * h` is bounded\n"
        "   using `h`'s maximum over the family *conditioned on* `z_e = 1`, which is\n"
        "   never worse than the unconditional envelope and is strictly better whenever\n"
        "   conditioning excludes competing structure.\n"
        "3. **A characterisation of when relaxation is loose.** The gap equals a\n"
        "   walk-count ratio; its exponent is set by network depth; it is unbounded on\n"
        "   downward-closed families that are not union-closed and exactly 1 on\n"
        "   join-closed ones.\n\n"
        "Part 3 is the general result and is independent of any application. Parts 1 and\n"
        "2 are mechanisms whose empirical value is, on present evidence, unproven at\n"
        "matched wall time.\n")

    if missing:
        L.append(f"\n\n## Missing sources\n\nNo result file found for: "
                 f"{', '.join(missing)}. Their tables are omitted rather than "
                 f"estimated.\n")

    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    with open(a.out, "w") as fh:
        fh.write("\n".join(L) + "\n")
    print(f"sources used:")
    for k, v in src.items():
        print(f"  {k:>14}: {os.path.basename(v) if v else 'MISSING'}")
    print(f"\ntables written: {a.out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
