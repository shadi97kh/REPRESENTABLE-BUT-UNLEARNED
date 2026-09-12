"""Exhaustive soundness audit of the affine bound propagation.

For every (family, model) in the config this enumerates the ENTIRE feasible family
independently and checks, at every feasible point:

  * every intermediate bracket   lo(z) <= true(z) <= hi(z)  at M1, H1, Z1, M2, H2, S, out
  * the numpy reference forward against the torch model
  * the final bound              U >= max_z f(X, z)
  * the incumbent                a real model evaluation at a feasible z, and <= U
  * combined vs terminal-only    U_combined <= U_terminal_only

A single violation aborts the run with a non-zero exit code and the offending case is
written out, because an unsound intermediate invalidates everything downstream.

Run:
    python -m scmp.audit_bounds --config configs/exhaustive.yaml
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from collections import Counter

import numpy as np
import torch
import yaml

from .bounds import (FLOAT_UNVERIFIED, INFEASIBLE, ResidualBounder,
                     certified_upper_bound)
from .model import CertifiedModel, backbone_adjacency
from .oracles import RNANonCrossingOracle
from .oracles.bruteforce import enumerate_rna


class Violation(Exception):
    def __init__(self, payload):
        super().__init__(json.dumps(payload, indent=2, default=str))
        self.payload = payload


def onehot(seq):
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return X


def build_model(seq, cands, d_hid, gamma, init_seed, scale, pseed):
    m = CertifiedModel(4, cands, len(seq), d_hid=d_hid, gamma=gamma, seed=init_seed)
    if scale:
        g = torch.Generator().manual_seed(pseed)
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype) * scale)
    return m


def zvec(m, s):
    z = np.zeros(len(m.candidates))
    for e in s:
        z[m.index[e]] = 1.0
    return z


def audit_case(cfg, fam, m, o, X, B, members, stats):
    slack = float(cfg["audit"]["slack"])
    seq = fam["seq"]

    if m.residual_active:
        rb = ResidualBounder(m, X, B)
        lo_aff, hi_aff, tr = rb.propagate()
        stats["envelope_inequalities"] += tr.envelope_terms
        for rows in tr.relu_kinds.values():
            for row in rows:
                stats["relu"].update(row)

        for s in members:
            z = zvec(m, s)
            truth = rb.reference_forward(z)

            if cfg["audit"].get("check_reference_forward", True):
                with torch.no_grad():
                    want = float(m.g(X, m.adjacency(torch.tensor(z, dtype=torch.float64), B)))
                got = float(truth["out"][0][0])
                stats["reference_checks"] += 1
                if abs(got - want) > 1e-9:
                    raise Violation({"kind": "reference_forward_mismatch", "seq": seq,
                                     "structure": s, "numpy": got, "torch": want})

            if cfg["audit"].get("check_intermediates", True):
                for name, (L, H) in tr.layers.items():
                    arr = truth[name]
                    for i in range(len(L)):
                        for k in range(len(L[i])):
                            lv, hv = L[i][k].eval(z), H[i][k].eval(z)
                            tv = float(arr[i][k])
                            stats["intermediate_checks"] += 1
                            if lv > tv + slack or tv > hv + slack:
                                raise Violation({
                                    "kind": "intermediate_bracket_violated",
                                    "seq": seq, "layer": name, "node": i, "unit": k,
                                    "structure": s, "lo": lv, "true": tv, "hi": hv})

    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)

    results = {}
    for mode in ("combined", "terminal_only"):
        r = certified_upper_bound(m, X, o, B, mode=mode)
        results[mode] = r
        stats["bounds"] += 1
        if r.upper < true_max - slack:
            raise Violation({"kind": "upper_bound_below_true_max", "seq": seq,
                             "mode": mode, "upper": r.upper, "true_max": true_max})
        with torch.no_grad():
            direct = float(m(X, m.indicator(r.incumbent_witness), B))
        if abs(direct - r.incumbent_value) > 0:
            raise Violation({"kind": "incumbent_is_not_a_model_evaluation", "seq": seq,
                             "mode": mode, "reported": r.incumbent_value,
                             "direct": direct})
        if tuple(sorted(r.incumbent_witness)) not in {tuple(sorted(x)) for x in members}:
            raise Violation({"kind": "incumbent_infeasible", "seq": seq, "mode": mode,
                             "witness": r.incumbent_witness})
        if r.incumbent_value > r.upper + slack:
            raise Violation({"kind": "incumbent_above_upper", "seq": seq, "mode": mode})

    rc, rt = results["combined"], results["terminal_only"]
    if rc.upper > rt.upper + slack:
        raise Violation({"kind": "combined_looser_than_separate_maxima", "seq": seq,
                         "combined": rc.upper, "terminal_only": rt.upper})
    if rc.model_hash != rt.model_hash:
        raise Violation({"kind": "baseline_weights_differ", "seq": seq})

    stats["strictly_tighter"] += rc.upper < rt.upper - 1e-9
    stats["gap_combined"].append(rc.upper - true_max)
    stats["gap_terminal"].append(rt.upper - true_max)
    if m.gamma == 0.0:
        stats["gamma0_cases"] += 1
        if abs(rc.upper - true_max) > 1e-12:
            raise Violation({"kind": "gamma_zero_not_exact", "seq": seq,
                             "upper": rc.upper, "true_max": true_max})
    return results


def audit_masks(cfg, m, o, X, B, seq, stats):
    mk = cfg.get("masks", {})
    g = o.ground_set()
    if mk.get("crossing_forced") and len(g) > 1:
        cross = [(a, b) for a in g for b in g if a[0] < b[0] < a[1] < b[1]]
        if cross:
            r = certified_upper_bound(m, X, o, B, forced=cross[0])
            stats["mask_checks"] += 1
            if r.termination_status != INFEASIBLE:
                raise Violation({"kind": "crossing_forced_not_infeasible", "seq": seq,
                                 "forced": cross[0]})
    if mk.get("forced_equals_forbidden") and g:
        r = certified_upper_bound(m, X, o, B, forced=(g[0],), forbidden=(g[0],))
        stats["mask_checks"] += 1
        if r.termination_status != INFEASIBLE:
            raise Violation({"kind": "contradictory_mask_not_infeasible", "seq": seq})
    slack = float(cfg["audit"]["slack"])
    for key, build in (("single_forced", lambda e: ((e,), ())),
                       ("single_forbidden", lambda e: ((), (e,)))):
        for e in g[:int(mk.get(key, 0))]:
            forced, forbidden = build(e)
            members = enumerate_rna(seq, 3, True, forced, forbidden)
            if not members:
                continue
            r = certified_upper_bound(m, X, o, B, forced=forced, forbidden=forbidden)
            stats["mask_checks"] += 1
            with torch.no_grad():
                tm = max(float(m(X, m.indicator(s), B)) for s in members)
            if r.upper < tm - slack:
                raise Violation({"kind": "masked_bound_below_true_max", "seq": seq,
                                 "forced": forced, "forbidden": forbidden,
                                 "upper": r.upper, "true_max": tm})


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.audit_bounds")
    ap.add_argument("--config", default="configs/exhaustive.yaml")
    args = ap.parse_args(argv)
    cfg = yaml.safe_load(open(args.config))

    t0 = time.time()
    stats = {"intermediate_checks": 0, "reference_checks": 0, "bounds": 0,
             "mask_checks": 0, "envelope_inequalities": 0, "strictly_tighter": 0,
             "gamma0_cases": 0, "cases": 0, "structures": 0,
             "relu": Counter(), "gap_combined": [], "gap_terminal": []}

    print(f"bound audit: {cfg['audit']['name']}  config={args.config}")
    print(f"slack={cfg['audit']['slack']:g}  stop_on_violation="
          f"{cfg['audit']['stop_on_violation']}\n")

    grid = list(itertools.product(cfg["models"]["d_hid"], cfg["models"]["gamma"],
                                 cfg["models"]["init_seed"],
                                 cfg["models"]["perturb_scale"],
                                 cfg["models"]["perturb_seed"]))
    try:
        for fam in cfg["families"]:
            seq = fam["seq"]
            o = RNANonCrossingOracle(seq, fam["min_loop"], fam["canonical_only"])
            members = enumerate_rna(seq, fam["min_loop"], fam["canonical_only"])
            X, B = onehot(seq), backbone_adjacency(len(seq))
            print(f"  {seq:14s} |E|={len(o.ground_set()):3d} |F|={len(members):5d} "
                  f"models={len(grid)}")
            for d_hid, gamma, iseed, scale, pseed in grid:
                m = build_model(seq, o.ground_set(), d_hid, gamma, iseed, scale, pseed)
                audit_case(cfg, fam, m, o, X, B, members, stats)
                audit_masks(cfg, m, o, X, B, seq, stats)
                stats["cases"] += 1
                stats["structures"] += len(members)
    except Violation as v:
        os.makedirs("runs/bound_violations", exist_ok=True)
        path = f"runs/bound_violations/{time.strftime('%Y%m%d-%H%M%S')}.json"
        with open(path, "w") as fh:
            json.dump(v.payload, fh, indent=2, default=str)
        print(f"\nSOUNDNESS VIOLATION -- aborting. Evidence: {path}")
        print(v)
        return 2

    gc = np.array(stats["gap_combined"])
    gt = np.array(stats["gap_terminal"])
    print(f"\ncases                     : {stats['cases']} "
          f"({stats['structures']} feasible structures visited)")
    print(f"intermediate bracket checks: {stats['intermediate_checks']}")
    print(f"reference-forward checks   : {stats['reference_checks']}")
    print(f"final bounds               : {stats['bounds']}")
    print(f"masked-bound checks        : {stats['mask_checks']}")
    print(f"envelope inequalities      : {stats['envelope_inequalities']} "
          f"(4 per product term)")
    print(f"relu units                 : {dict(stats['relu'])}")
    print(f"gamma=0 cases exact        : {stats['gamma0_cases']}")
    print(f"combined strictly tighter  : {stats['strictly_tighter']}/{len(gc)}")
    if len(gc):
        print(f"gap to true max, combined  : median {np.median(gc):.4f}  "
              f"max {gc.max():.4f}")
        print(f"gap to true max, terminal  : median {np.median(gt):.4f}  "
              f"max {gt.max():.4f}")
    print(f"\nnumerical status           : {FLOAT_UNVERIFIED} -- DIAGNOSTIC, not a proof")
    print(f"                             see docs/bounds_proof.md")
    print(f"runtime                    : {time.time() - t0:.2f}s")
    print("\nno soundness violation found")
    return 0


if __name__ == "__main__":
    sys.exit(main())
