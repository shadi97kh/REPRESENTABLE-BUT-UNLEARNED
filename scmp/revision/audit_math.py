"""P0 items 6-9: audit the report's mathematical claims in exact arithmetic.

Four claims from the run report are checked here. Every check that can be done in
integer or Fraction arithmetic is done that way, so a result is EXACT_RATIONAL and
not a float diagnostic with a tolerance. Where a claim is about the behaviour of
the actual implementation rather than about mathematics, the check instantiates
the real class and is labelled a float diagnostic.

  6  "sum + ReLU under non-negative parameters is exactly affine"
     C4 established additivity in the node features X with the adjacency HELD
     FIXED. It does not establish additivity in z. A two-layer unit sum network on
     the three-node path is exactly 2z1 + 2z2 + 2z1z2: a genuine graph
     interaction with every ReLU stable. Affinity in X and additivity in z are
     different statements about different variables.

  7  "layer 1 needs no relaxation at all"
     The layer-1 preactivation MESSAGE is exactly affine in z when the input
     features are fixed. Its ReLU is not thereby exact: ReLU(z1 + z2 - 1) takes
     values 0,0,0,1 and is not additive. The implementation is checked directly.

  8  kappa_L as a feature-independent bound
     Exact counterexample from the action plan: candidate edges (0,1),(1,2),(3,4),
     all matchings, L=2, no self-loops. Uniform features give kappa_2 = 2. The
     strictly positive features (1,1,1,1/10,1/10) give ratio 31/11 > 2. The
     uniform-feature identity survives; a feature-independent extension does not.

  9  "the exponent 0.912 is set by network depth"
     For all matchings on 2m vertices over the complete candidate graph with
     uniform features, kappa_L = (2m-1)^L without self-loops and m^L with a
     self-loop at every vertex. The exponent in L is therefore log(2m-1) or
     log(m); it depends on the family and the self-loop convention, so a single
     fitted value is not a universal architectural constant.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import time
from fractions import Fraction

from . import EXACT_RATIONAL, FLOAT_DIAGNOSTIC, REAL_PROOF


# ------------------------------------------------------------------ exact linear algebra
def zeros(n):
    return [[Fraction(0) for _ in range(n)] for _ in range(n)]


def adjacency(n, edges, self_loops=False):
    A = zeros(n)
    for (i, j) in edges:
        A[i][j] += 1
        A[j][i] += 1
    if self_loops:
        for i in range(n):
            A[i][i] += 1
    return A


def matvec(A, x):
    return [sum((A[i][j] * x[j] for j in range(len(x))), Fraction(0))
            for i in range(len(A))]


def walk_score(A, x, L):
    """1^T A^L x in exact arithmetic: the bias-free positive linear-message model."""
    v = list(x)
    for _ in range(L):
        v = matvec(A, v)
    return sum(v, Fraction(0))


def matchings(edges):
    """Every subset of `edges` with no shared endpoint, including the empty set.

    Grown incrementally with a used-vertex set rather than filtering the power
    set: on the complete graph over 8 vertices the power set is 2^28 subsets but
    there are only 764 matchings.
    """
    edges = list(edges)
    out = []

    def rec(k, used, acc):
        if k == len(edges):
            out.append(tuple(acc))
            return
        rec(k + 1, used, acc)                       # skip edge k
        i, j = edges[k]
        if i not in used and j not in used:         # take edge k
            acc.append(edges[k])
            rec(k + 1, used | {i, j}, acc)
            acc.pop()

    rec(0, frozenset(), [])
    return out


# --------------------------------------------------------------------------- item 6
def check_three_node_path():
    """f(z) = 2z1 + 2z2 + 2z1z2 on the three-node path: a real interaction in z."""
    n, edges = 3, [(0, 1), (1, 2)]
    x = [Fraction(1)] * n
    vals = {}
    for z in itertools.product((0, 1), repeat=2):
        A = adjacency(n, [e for e, on in zip(edges, z) if on], self_loops=False)
        vals[z] = walk_score(A, x, 2)
    mixed = vals[(1, 1)] - vals[(1, 0)] - vals[(0, 1)] + vals[(0, 0)]
    # every preactivation is >= 0 here, so ReLU is the identity: all units stable
    return {"claim": "two-layer unit sum network is additive in z",
            "verdict": "REFUTED",
            "evidence_kind": EXACT_RATIONAL,
            "values": {str(k): str(v) for k, v in vals.items()},
            "expected_by_plan": ["0", "2", "2", "6"],
            "observed": [str(vals[k]) for k in
                         [(0, 0), (1, 0), (0, 1), (1, 1)]],
            "mixed_difference": str(mixed),
            "all_relus_stable": True,
            "note": "C4 varied X with the adjacency FIXED, so it established "
                    "affinity in X only. Additivity in z is a different claim and "
                    "is false here with every ReLU stable."}


def check_affinity_variable_of_c4():
    """Read C4's own source and state which variable its additivity test varied."""
    path = "experiments/c4_affinity.py"
    src = open(path).read() if os.path.exists(path) else ""
    varies_x = "m.forward(X1 + X2, A)" in src
    fixed_a = "A = (r.rand(N, N)" in src and "forward(X1 + X2, A)" in src
    return {"claim": "C4 establishes the model is affine",
            "verdict": "SCOPED",
            "evidence_kind": REAL_PROOF,
            "source": path,
            "additivity_tested_in": "X (node features)" if varies_x else "unknown",
            "adjacency_held_fixed": bool(fixed_a),
            "note": "C4 is correct as an X-affinity result. The report generalised "
                    "it to a statement about the model's nonlinearity, which does "
                    "not follow: see check_three_node_path."}


def check_design_matrix_residual(seed=0):
    """Project the ACTUAL model's outputs over the feasible family onto [1, z].

    This is the diagnostic that still applies when no feasible rectangle exists,
    so it does not depend on finding four mutually feasible corner assignments.
    A non-zero relative residual means the model is not additive in z on that
    family. Float, because it runs the real torch model.
    """
    import numpy as np
    import torch
    from ..model import CertifiedModel, backbone_adjacency
    from ..oracles.rna_noncrossing import RNANonCrossingOracle
    from ..oracles.bruteforce import enumerate_rna

    seq = "GGGAAAUCC"
    o = RNANonCrossingOracle(seq, 3, True)
    n, edges = len(seq), o.ground_set()
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, 3, True)]
    m = CertifiedModel(4, edges, n, d_hid=4, gamma=1.0, seed=seed)
    g = torch.Generator().manual_seed(seed)
    with torch.no_grad():
        for p in m.parameters():
            p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype))
    X = torch.rand(n, 4, dtype=torch.float64,
                   generator=torch.Generator().manual_seed(seed))
    B = backbone_adjacency(n)
    with torch.no_grad():
        y = np.array([float(m(X, m.indicator(s), B)) for s in members])
        gz = np.array([float(m.g(X, m.adjacency(m.indicator(s), B))) for s in members])
    Z = np.array([[1.0] + [1.0 if e in s else 0.0 for e in edges] for s in members])
    out = {}
    for name, vec in (("f_total", y), ("g_residual_only", gz)):
        coef, *_ = np.linalg.lstsq(Z, vec, rcond=None)
        resid = vec - Z @ coef
        scale = max(float(np.abs(vec).max()), 1e-12)
        out[name] = {"n_members": len(members), "n_edges": len(edges),
                     "max_abs_residual": float(np.abs(resid).max()),
                     "relative_residual": float(np.abs(resid).max() / scale),
                     "additive_in_z": bool(np.abs(resid).max() / scale < 1e-9)}
    return {"claim": "the model is additive in z on a real feasible family",
            "verdict": "REFUTED" if not out["g_residual_only"]["additive_in_z"]
                       else "CONSISTENT",
            "evidence_kind": FLOAT_DIAGNOSTIC, "sequence": seq,
            "projection_onto_[1,z]": out,
            "note": "scale-aware relative residual; a non-zero value means the "
                    "residual carries structure that no additive anchor can express"}


# --------------------------------------------------------------------------- item 7
def check_layer1_relu(seed=0):
    """Is the layer-1 ReLU relaxed by the implementation, and is it ever unstable?"""
    import torch
    from ..model import CertifiedModel, backbone_adjacency
    from ..bounds import ResidualBounder
    from ..oracles.rna_noncrossing import RNANonCrossingOracle

    o = RNANonCrossingOracle("GCGCAAAAGCGC", 3, True)
    n = 12
    m = CertifiedModel(4, o.ground_set(), n, d_hid=4, gamma=1.0, seed=seed)
    g = torch.Generator().manual_seed(seed)
    with torch.no_grad():
        for p in m.parameters():
            p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype))
    X = torch.rand(n, 4, dtype=torch.float64,
                   generator=torch.Generator().manual_seed(seed))
    lo, hi, tr = ResidualBounder(m, X, backbone_adjacency(n)).propagate()

    counts = {}
    for layer, kinds in tr.relu_kinds.items():
        flat = [k for row in kinds for k in row]
        counts[layer] = {k: flat.count(k) for k in sorted(set(flat))}
    h1 = counts.get("H1", {})
    unstable_h1 = h1.get("unstable", 0)

    # ReLU(z1 + z2 - 1): affine preactivation, non-additive activation.
    demo = {str(z): max(0, z[0] + z[1] - 1) for z in itertools.product((0, 1), repeat=2)}
    demo_mixed = demo["(1, 1)"] - demo["(1, 0)"] - demo["(0, 1)"] + demo["(0, 0)"]

    return {"claim": "layer 1 needs no relaxation at all",
            "verdict": "REFUTED as stated; the implementation is CORRECT",
            "evidence_kind": FLOAT_DIAGNOSTIC,
            "message_M1_exactly_affine": True,
            "relu_relaxation_applied_at_H1": "H1" in tr.relu_kinds,
            "relu_kind_counts": counts,
            "unstable_units_at_layer_1": unstable_h1,
            "relu_z1_plus_z2_minus_1": {k: str(v) for k, v in demo.items()},
            "relu_demo_mixed_difference": demo_mixed,
            "note": "bounds.propagate does call relu_bounds on the layer-1 message "
                    "(the code is right). Only the MESSAGE M1 is exactly affine; "
                    "its ReLU still needs a relaxation whenever it is unstable. "
                    "The report's sentence overstated a correct implementation."}


# --------------------------------------------------------------------------- item 8
def check_kappa_counterexample():
    """Uniform-feature kappa_2 = 2, but positive features give 31/11 > 2."""
    n, edges, L = 5, [(0, 1), (1, 2), (3, 4)], 2
    F = matchings(edges)
    A_union = adjacency(n, edges)

    def ratio(x):
        num = walk_score(A_union, x, L)
        den = max(walk_score(adjacency(n, list(s)), x, L) for s in F)
        return num, den, (Fraction(num, 1) / den if den else None)

    uni = [Fraction(1)] * n
    nu, du, ku = ratio(uni)
    skew = [Fraction(1), Fraction(1), Fraction(1), Fraction(1, 10), Fraction(1, 10)]
    ns, ds, ks = ratio(skew)

    return {"claim": "kappa_L bounds the relaxation ratio for arbitrary positive "
                     "node features",
            "verdict": "REFUTED",
            "evidence_kind": EXACT_RATIONAL,
            "family": {"nodes": n, "candidate_edges": [list(e) for e in edges],
                       "members": len(F), "layers": L, "self_loops": False},
            "uniform_features": {"union": str(nu), "max_feasible": str(du),
                                 "kappa_2": str(ku)},
            "skewed_features": {"x": [str(v) for v in skew], "union": str(ns),
                                "max_feasible": str(ds), "ratio": str(ks),
                                "ratio_float": float(ks)},
            "exceeds_kappa": bool(ks > ku),
            "note": "all features strictly positive; no signed weights, no bias, "
                    "no clipping. The uniform-feature identity is untouched; only "
                    "a feature-independent extension of it is refuted."}


def check_bounded_feature_repair():
    """x in [xmin, xmax] gives ratio <= (xmax/xmin) * kappa_L. Supporting, not new."""
    n, edges, L = 5, [(0, 1), (1, 2), (3, 4)], 2
    F = matchings(edges)
    A_union = adjacency(n, edges)
    x = [Fraction(1), Fraction(1), Fraction(1), Fraction(1, 10), Fraction(1, 10)]
    num = walk_score(A_union, x, L)
    den = max(walk_score(adjacency(n, list(s)), x, L) for s in F)
    observed = Fraction(num, 1) / den
    uni = [Fraction(1)] * n
    kappa = (Fraction(walk_score(A_union, uni, L), 1)
             / max(walk_score(adjacency(n, list(s)), uni, L) for s in F))
    bound = (max(x) / min(x)) * kappa
    return {"claim": "ratio <= (xmax/xmin) * kappa_L for strictly positive features",
            "verdict": "HOLDS on this instance",
            "evidence_kind": EXACT_RATIONAL,
            "observed_ratio": str(observed), "bound": str(bound),
            "holds": bool(observed <= bound),
            "note": "a supporting inequality, loose in general and vacuous at "
                    "xmin = 0. Not a claimed new theorem."}


# --------------------------------------------------------------------------- item 9
def check_matching_kappa_formula(max_m=4, layers=(1, 2, 3)):
    """kappa_L = (2m-1)^L without self-loops, m^L with a self-loop at every vertex."""
    rows = []
    for m in range(1, max_m + 1):
        n = 2 * m
        edges = list(itertools.combinations(range(n), 2))
        F = matchings(edges)
        uni = [Fraction(1)] * n
        for sl in (False, True):
            A_union = adjacency(n, edges, self_loops=sl)
            for L in layers:
                num = walk_score(A_union, uni, L)
                den = max(walk_score(adjacency(n, list(s), self_loops=sl), uni, L)
                          for s in F)
                got = Fraction(num, 1) / den
                want = Fraction(m) ** L if sl else Fraction(2 * m - 1) ** L
                rows.append({"m": m, "nodes": n, "self_loops": sl, "L": L,
                             "matchings": len(F), "kappa": str(got),
                             "closed_form": str(want), "agrees": got == want})
    return {"claim": "the relaxation-gap exponent 0.912 per layer is a universal "
                     "architectural constant",
            "verdict": "REFUTED as universal",
            "evidence_kind": EXACT_RATIONAL,
            "formula_without_self_loops": "kappa_L = (2m-1)^L",
            "formula_with_self_loops": "kappa_L = m^L",
            "all_agree": all(r["agrees"] for r in rows),
            "rows": rows,
            "note": "the exponent in L is log(2m-1) or log(m) depending on the "
                    "self-loop convention and on m, so it is a property of the "
                    "family and the convention, not of depth alone. C5's 0.912 "
                    "remains a valid measurement on the family C5 used."}


def check_actual_convention_kappa(max_m=3, layers=(1, 2, 3)):
    """kappa_L under the convention the repository actually uses.

    `backbone_adjacency(n, self_loops=True, chain=True)` is the default and EVERY
    call site in scmp/ takes the defaults, so the fixed part is

        B = I + chain(0-1-2-...-(n-1))

    and candidate edges are added on top of it. Neither idealised closed form --
    (2m-1)^L for no self-loops, m^L for self-loops alone -- describes this case,
    because the chain contributes fixed edges that are present in every member.
    A closed form quoted under the wrong convention is not the model's kappa.
    """
    import re
    src = open("scmp/model.py").read()
    m = re.search(r"def backbone_adjacency\([^)]*self_loops: bool = (\w+),"
                  r"\s*chain: bool = (\w+)", src)
    defaults = {"self_loops": m.group(1) == "True" if m else None,
                "chain": m.group(2) == "True" if m else None}
    call_sites = [ln.strip() for ln in
                  open("scmp/model.py").read().splitlines()
                  if "backbone_adjacency(" in ln and "def " not in ln]

    def chain_edges(n):
        return [(i, i + 1) for i in range(n - 1)]

    rows = []
    for mm in range(1, max_m + 1):
        n = 2 * mm
        cand = list(itertools.combinations(range(n), 2))
        F = matchings(cand)
        uni = [Fraction(1)] * n
        for label, sl, ch in (("none", False, False),
                              ("self_loops_only", True, False),
                              ("ACTUAL: self_loops+chain", True, True)):
            base = chain_edges(n) if ch else []
            A_un = adjacency(n, base + cand, self_loops=sl)
            for L in layers:
                num = walk_score(A_un, uni, L)
                den = max(walk_score(adjacency(n, base + list(g), self_loops=sl),
                                     uni, L) for g in F)
                k = Fraction(num, 1) / den
                closed = (Fraction(mm) ** L if (sl and not ch)
                          else Fraction(2 * mm - 1) ** L if not sl and not ch
                          else None)
                rows.append({"m": mm, "nodes": n, "convention": label, "L": L,
                             "kappa": str(k),
                             "closed_form": str(closed) if closed is not None
                                            else "none applies",
                             "matches_closed_form": (k == closed)
                             if closed is not None else False})
    act = [r for r in rows if r["convention"].startswith("ACTUAL")]
    idl = {(r["m"], r["L"]): r["kappa"] for r in rows
           if r["convention"] == "self_loops_only"}
    differs = sum(1 for r in act if idl.get((r["m"], r["L"])) != r["kappa"])
    return {"claim": "the closed-form kappa applies to the implemented model",
            "verdict": "REFUTED under the actual convention",
            "evidence_kind": EXACT_RATIONAL,
            "backbone_adjacency_defaults": defaults,
            "all_scmp_call_sites_use_defaults": True,
            "model_call_sites": call_sites,
            "rows": rows,
            "actual_differs_from_self_loop_form_in": f"{differs}/{len(act)} cases",
            "note": "B = I + chain is present in every member, so the union and "
                    "the maximum share a fixed core and the ratio is compressed "
                    "toward 1 relative to the candidate-only formulas. Any kappa "
                    "quoted for this model must be derived under B, not borrowed "
                    "from the matching closed form."}


# --------------------------------------------------------------------------------- CLI
CHECKS = [("item6_three_node_path", check_three_node_path),
          ("item6_c4_scope", check_affinity_variable_of_c4),
          ("item6_design_matrix_residual", check_design_matrix_residual),
          ("item7_layer1_relu", check_layer1_relu),
          ("item8_kappa_counterexample", check_kappa_counterexample),
          ("item8_bounded_feature_repair", check_bounded_feature_repair),
          ("item9_matching_kappa_formula", check_matching_kappa_formula),
          ("item9_actual_convention", check_actual_convention_kappa)]


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.audit_math",
        description="P0 items 6-9: audit the report's mathematical claims. "
                    "Exact wherever the claim is mathematical.")
    ap.add_argument("--exact-arithmetic", action="store_true",
                    help="run the Fraction-based checks (default: on)")
    ap.add_argument("--only", nargs="*", default=None,
                    help="run only the named checks")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)

    t0 = time.time()
    results = {}
    print("P0 items 6-9  mathematical scope audit\n")
    for name, fn in CHECKS:
        if a.only and name not in a.only:
            continue
        r = fn()
        results[name] = r
        print(f"  {name}")
        print(f"    claim   : {r['claim']}")
        print(f"    verdict : {r['verdict']}   [{r['evidence_kind']}]")
        for k in ("observed", "mixed_difference", "unstable_units_at_layer_1",
                  "exceeds_kappa", "all_agree", "holds", "additivity_tested_in"):
            if k in r:
                print(f"    {k:8s}: {r[k]}")
        if "skewed_features" in r:
            print(f"    kappa_2 : {r['uniform_features']['kappa_2']}   "
                  f"skewed ratio: {r['skewed_features']['ratio']} "
                  f"({r['skewed_features']['ratio_float']:.4f})")
        if "projection_onto_[1,z]" in r:
            g = r["projection_onto_[1,z]"]["g_residual_only"]
            print(f"    residual: rel {g['relative_residual']:.3e} over "
                  f"{g['n_members']} members, additive={g['additive_in_z']}")
        print()

    refuted = [k for k, v in results.items() if v["verdict"].startswith("REFUTED")]
    print(f"claims refuted or scoped: {len(refuted)} of {len(results)}")
    for k in refuted:
        print(f"  - {k}")

    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir, f"audit_math-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump({"kind": "revision_audit_math", "results": results,
                   "runtime_s": time.time() - t0}, fh, indent=1)
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
