"""Model-effect intervals for a sub-family D contained in F.

THE INTERVAL. With sound bounds L_F <= max_F f <= U_F and L_D <= max_D f <= U_D, and
D contained in F (so max_D f <= max_F f, hence the effect is non-negative):

    effect := max_F f - max_D f   lies in   [ max(0, L_F - U_D),  U_F - L_D ]

Both endpoints are needed. The lower endpoint uses the LOWER bound on F against the UPPER
bound on D; the upper endpoint the other way round. Each L is an incumbent, i.e. a value
the model actually attains, so both endpoints are certified rather than estimated.

WHY U_F - U_D IS NOT THE EFFECT. Both are upper bounds carrying unknown and different
slack. Their difference can be positive when the true effect is zero, or negative when it
is large. It is not a bound in either direction and this module refuses to report it as
one; it is printed only in a column explicitly labelled as not an effect.

SEPARATION OF CLAIMS. Three things are kept apart and never merged:
  * model-effect intervals   -- statements about the model's output over a family;
  * threshold decisions      -- whether a certified bound clears a stated tau;
  * biological causal claims -- none are made here, and none are supportable: the
                                label is marginal over latent structures
                                (preregistration.md, section 6).

ABSTENTION. When either side is unresolved, or the interval contains zero and is wider
than the decision margin, the verdict is `abstain`. An abstention is a result.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys
import time
from dataclasses import dataclass

import numpy as np
import torch
import yaml

from .bounds import certified_upper_bound
from .model import CertifiedModel, backbone_adjacency
from .oracles import PathMatchingOracle, RNANonCrossingOracle
from .oracles.bruteforce import enumerate_rna


@dataclass
class ModelEffect:
    edge: tuple
    L_F: float
    U_F: float
    L_D: float
    U_D: float
    lo: float
    hi: float
    verdict: str
    width: float
    note: str = ""

    @property
    def not_an_effect_UF_minus_UD(self) -> float:
        """Printed for contrast only. Never a bound on the effect."""
        return self.U_F - self.U_D


def model_effect_interval(L_F, U_F, L_D, U_D, margin=1e-6) -> ModelEffect:
    lo = max(0.0, L_F - U_D)
    hi = U_F - L_D
    if not (np.isfinite(lo) and np.isfinite(hi)):
        verdict = "abstain"
    elif hi < lo:
        verdict = "invalid"          # would mean an unsound input bound
    elif lo > margin:
        verdict = "effect_positive"
    elif hi <= margin:
        verdict = "effect_negligible"
    else:
        verdict = "abstain"
    return ModelEffect((), L_F, U_F, L_D, U_D, lo, hi, verdict, hi - lo)


def effects_for_instance(m, X, o, B, edges, margin):
    """One interval per conditioning edge e, with D = F and {z_e = 1}."""
    rF = certified_upper_bound(m, X, o, B, mode="combined")
    out = []
    for e in edges:
        feas = o.feasibility(forced=(e,))
        if not feas.feasible:
            out.append(ModelEffect(e, rF.incumbent_value, rF.upper, float("-inf"),
                                   float("-inf"), 0.0, 0.0, "subfamily_empty", 0.0,
                                   feas.reason))
            continue
        rD = certified_upper_bound(m, X, o, B, forced=(e,), mode="combined")
        me = model_effect_interval(rF.incumbent_value, rF.upper,
                                   rD.incumbent_value, rD.upper, margin)
        out.append(ModelEffect(e, me.L_F, me.U_F, me.L_D, me.U_D, me.lo, me.hi,
                               me.verdict, me.width))
    return rF, out


def validate_by_enumeration(m, X, o, B, members, edges, margin):
    """On tiny families the true effect is computable; every interval must contain it."""
    with torch.no_grad():
        fvals = {tuple(sorted(s)): float(m(X, m.indicator(s), B)) for s in members}
    trueF = max(fvals.values())
    fails, checked = [], 0
    _, effects = effects_for_instance(m, X, o, B, edges, margin)
    for me in effects:
        if me.verdict == "subfamily_empty":
            continue
        sub = [v for s, v in fvals.items() if me.edge in s]
        if not sub:
            continue
        true_effect = trueF - max(sub)
        checked += 1
        if not (me.lo - 1e-9 <= true_effect <= me.hi + 1e-9):
            fails.append({"edge": list(me.edge), "true_effect": true_effect,
                          "interval": [me.lo, me.hi]})
        # the contrast quantity: show it is NOT a bound
    return checked, fails, effects, trueF


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.evaluate_explanations")
    ap.add_argument("--config", default="configs/pilot.yaml")
    ap.add_argument("--margin", type=float, default=1e-6)
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    t0 = time.time()

    print("model-effect intervals   [max(0, L_F - U_D),  U_F - L_D]")
    print("U_F - U_D is NOT an effect and is shown only for contrast\n")

    results = {"validation": [], "instances": [], "contrast_violations": 0}

    # ---- validation on tiny families, by enumeration ----------------------
    print("validation by exhaustive enumeration on tiny families")
    total_checked, total_fails = 0, []
    cases = [("path_matching", PathMatchingOracle(9), 9),
             ("rna_noncrossing", RNANonCrossingOracle("GGGAAAUCC", 3, True), 9),
             ("rna_noncrossing", RNANonCrossingOracle("GCGCAAAAGCGC", 3, True), 12)]
    for name, o, n in cases:
        edges = o.ground_set()
        if name == "path_matching":
            members = [s for r in range(len(edges) + 1)
                       for s in itertools.combinations(edges, r)
                       if all(not (set(x) & set(y))
                              for i, x in enumerate(s) for y in s[i + 1:])]
        else:
            members = enumerate_rna(o.seq, o.min_loop, o.canonical_only)
        for seed in (0, 1, 2):
            m = CertifiedModel(4, edges, n, d_hid=4, gamma=1.0, seed=seed)
            g = torch.Generator().manual_seed(100 + seed)
            with torch.no_grad():
                for p in m.parameters():
                    p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype))
            X = torch.rand(n, 4, dtype=torch.float64,
                           generator=torch.Generator().manual_seed(seed))
            B = backbone_adjacency(n)
            checked, fails, effects, trueF = validate_by_enumeration(
                m, X, o, B, members, edges, a.margin)
            total_checked += checked
            total_fails += fails
            # count how often U_F - U_D would have been a wrong "effect"
            with torch.no_grad():
                fvals = {tuple(sorted(s)): float(m(X, m.indicator(s), B))
                         for s in members}
            for me in effects:
                if me.verdict == "subfamily_empty":
                    continue
                sub = [v for s, v in fvals.items() if me.edge in s]
                if not sub:
                    continue
                te = max(fvals.values()) - max(sub)
                naive = me.not_an_effect_UF_minus_UD
                if abs(naive - te) > max(me.width, 1e-9):
                    results["contrast_violations"] += 1
        print(f"  {name:16s} n={n:3d} |F|={len(members):5d}  intervals checked "
              f"{checked * 3:5d}")
    print(f"\n  intervals validated : {total_checked}")
    print(f"  containment failures: {len(total_fails)}")
    if total_fails:
        for f in total_fails[:5]:
            print(f"    {f}")
        print("\nSOUNDNESS VIOLATION in the model-effect interval; aborting")
        return 2
    print(f"  cases where U_F - U_D would have misstated the effect by more than the "
          f"interval width: {results['contrast_violations']}")

    # ---- report on pilot-scale instances ---------------------------------
    print("\nmodel-effect intervals on pilot instances (margin "
          f"{a.margin:g})")
    verdicts = {}
    o = RNANonCrossingOracle("GCGCAAAAGCGC", 3, True)
    n = 12
    rows = []
    for seed in cfg["seeds"]:
        m = CertifiedModel(4, o.ground_set(), n, d_hid=4, gamma=1.0, seed=seed)
        g = torch.Generator().manual_seed(seed)
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype))
        X = torch.rand(n, 4, dtype=torch.float64,
                       generator=torch.Generator().manual_seed(seed))
        B = backbone_adjacency(n)
        rF, effects = effects_for_instance(m, X, o, B, o.ground_set(), a.margin)
        for me in effects:
            verdicts[me.verdict] = verdicts.get(me.verdict, 0) + 1
            rows.append({"seed": seed, "edge": list(me.edge), "lo": me.lo,
                         "hi": me.hi, "width": me.width, "verdict": me.verdict,
                         "L_F": me.L_F, "U_F": me.U_F, "L_D": me.L_D, "U_D": me.U_D,
                         "not_an_effect_UF_minus_UD": me.not_an_effect_UF_minus_UD})
    print(f"  {'verdict':>20} {'count':>6}")
    for k, v in sorted(verdicts.items(), key=lambda t: -t[1]):
        print(f"  {k:>20} {v:6d}")
    widths = [r["width"] for r in rows]
    print(f"\n  median interval width: {np.median(widths):.4f}")
    print(f"  abstention rate      : "
          f"{verdicts.get('abstain', 0) / max(len(rows), 1):.1%}")

    print("\nseparation of claims")
    print("  model-effect intervals : reported above, model-internal")
    print("  threshold decisions    : NOT computed here; a tau decision is a separate")
    print("                           artefact and is not implied by any interval")
    print("  paired interventions   : NOT computed here")
    print("  biological causal      : NONE. The label is marginal over latent")
    print("                           structures (preregistration.md section 6), so no")
    print("                           causal or efficacy claim follows from these")

    out = "runs/explanations"
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"effects-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump({"validated_intervals": total_checked,
                   "containment_failures": len(total_fails),
                   "contrast_violations": results["contrast_violations"],
                   "verdicts": verdicts, "rows": rows,
                   "margin": a.margin,
                   "numerical_status": "float_unverified",
                   "claims": {"biological_causal": "none",
                              "threshold_decisions": "separate artefact",
                              "paired_interventions": "not computed"}},
                  fh, indent=2)
    print(f"\nwritten: {path}")
    print(f"runtime: {time.time() - t0:.1f}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
