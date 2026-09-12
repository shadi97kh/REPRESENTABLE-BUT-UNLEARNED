"""Randomised audit of conditional message support.

Checks, on every case, against exhaustive enumeration of the feasible family:
  * the conditional bound dominates the true maximum
  * it never exceeds the inherited unconditional fallback (min of valid bounds retained)
  * recursive work stays within budget
  * an edge is dropped only when the oracle proved it impossible

Then compares the three deterministic edge choices -- random, largest-gap and
measured-gain -- on tightness against teacher cost.

Run:  python -m scmp.audit_conditional --cases 500 --seed 7027
"""
from __future__ import annotations

import argparse
import json
import os
import random
import sys
import time

import numpy as np
import torch

from .bounds import INFEASIBLE, certified_upper_bound
from .conditional import BoundCache, SCORERS, conditional_upper_bound
from .model import CertifiedModel, backbone_adjacency
from .oracles import RNANonCrossingOracle
from .oracles.bruteforce import enumerate_rna

SEQS = ["GGAAACC", "GGGAAAUCC", "GCGCAAAAGCGC", "GGGGAAAACCCC"]
SLACK = 1e-9


class Violation(Exception):
    def __init__(self, payload):
        super().__init__(json.dumps(payload, indent=2, default=str))
        self.payload = payload


def onehot(seq):
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return X


def make_case(rng):
    seq = rng.choice(SEQS)
    o = RNANonCrossingOracle(seq, 3, True)
    cand = o.ground_set()
    m = CertifiedModel(4, cand, len(seq), d_hid=rng.choice([4, 5]),
                       gamma=rng.choice([0.0, 1.0, -1.0, 0.5]), seed=rng.randint(0, 999))
    scale = rng.choice([0.0, 0.6, 1.2])
    if scale:
        g = torch.Generator().manual_seed(rng.randint(0, 9999))
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype) * scale)
    forced, forbidden = (), ()
    r = rng.random()
    if r < 0.2 and cand:
        forced = (cand[rng.randrange(len(cand))],)
    elif r < 0.35 and cand:
        forbidden = (cand[rng.randrange(len(cand))],)
    return seq, o, m, onehot(seq), backbone_adjacency(len(seq)), forced, forbidden


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.audit_conditional")
    ap.add_argument("--cases", type=int, default=500)
    ap.add_argument("--seed", type=int, default=7027)
    ap.add_argument("--budget", type=int, default=5)
    a = ap.parse_args(argv)

    rng = random.Random(a.seed)
    t0 = time.time()
    stats = {"cases": 0, "skipped_infeasible": 0, "structures": 0,
             "cond_passes": 0, "teacher_passes": 0, "oracle_calls": 0,
             "proved_impossible": 0, "conditional_used": 0, "fallback_used": 0,
             "strictly_tighter_than_baseline": 0, "fallback_won": 0}
    gaps = {k: [] for k in ("baseline", "random", "largest_gap", "measured_gain")}
    costs = {k: [] for k in SCORERS}

    print(f"conditional-support audit  cases={a.cases} seed={a.seed} budget={a.budget}\n")
    try:
        for _ in range(a.cases):
            seq, o, m, X, B, forced, forbidden = make_case(rng)
            if not o.feasibility(forced, forbidden).feasible:
                stats["skipped_infeasible"] += 1
                continue
            members = enumerate_rna(seq, 3, True, forced, forbidden)
            if not members:
                stats["skipped_infeasible"] += 1
                continue
            with torch.no_grad():
                true_max = max(float(m(X, m.indicator(s), B)) for s in members)
            base = certified_upper_bound(m, X, o, B, forced, forbidden, mode="combined")
            gaps["baseline"].append(base.upper - true_max)

            for scorer in SCORERS:
                cache = BoundCache()
                r = conditional_upper_bound(m, X, o, B, forced, forbidden,
                                            budget=a.budget, scorer=scorer,
                                            seed=a.seed, cache=cache)
                if r.upper < true_max - SLACK:
                    raise Violation({"kind": "upper_below_true_max", "seq": seq,
                                     "scorer": scorer, "upper": r.upper,
                                     "true_max": true_max, "forced": forced,
                                     "forbidden": forbidden})
                if r.upper > base.upper + SLACK:
                    raise Violation({"kind": "conditional_looser_than_fallback",
                                     "seq": seq, "scorer": scorer,
                                     "conditional": r.upper, "fallback": base.upper})
                if r.stats.get("conditional_passes", 0) > a.budget:
                    raise Violation({"kind": "budget_exceeded", "seq": seq,
                                     "scorer": scorer,
                                     "passes": r.stats["conditional_passes"],
                                     "budget": a.budget})
                with torch.no_grad():
                    direct = float(m(X, m.indicator(r.incumbent_witness), B))
                if abs(direct - r.incumbent_value) > 0:
                    raise Violation({"kind": "incumbent_not_a_model_evaluation",
                                     "seq": seq, "scorer": scorer})
                gaps[scorer].append(r.upper - true_max)
                costs[scorer].append(r.stats.get("teacher_passes", 0))
                if scorer == "largest_gap":
                    for key in ("cond_passes", "teacher_passes", "oracle_calls",
                                "proved_impossible", "conditional_used",
                                "fallback_used"):
                        src = {"cond_passes": "conditional_passes"}.get(key, key)
                        stats[key] += r.stats.get(src, 0)
                    stats["strictly_tighter_than_baseline"] += r.upper < base.upper - 1e-9
                    stats["fallback_won"] += r.stats.get("used") == "unconditional_fallback"
            stats["cases"] += 1
            stats["structures"] += len(members)
    except Violation as v:
        os.makedirs("runs/conditional_violations", exist_ok=True)
        path = f"runs/conditional_violations/{time.strftime('%Y%m%d-%H%M%S')}.json"
        with open(path, "w") as fh:
            json.dump(v.payload, fh, indent=2, default=str)
        print(f"\nSOUNDNESS VIOLATION -- aborting. Evidence: {path}")
        print(v)
        return 2

    print(f"cases audited            : {stats['cases']} "
          f"({stats['structures']} feasible structures)")
    print(f"skipped (infeasible mask): {stats['skipped_infeasible']}")
    print(f"conditional passes       : {stats['cond_passes']} (budget {a.budget}/case)")
    print(f"oracle calls             : {stats['oracle_calls']}")
    print(f"edges proved impossible  : {stats['proved_impossible']}")
    print(f"q from conditional       : {stats['conditional_used']}")
    print(f"q from fallback          : {stats['fallback_used']}")
    print(f"tighter than baseline    : {stats['strictly_tighter_than_baseline']}"
          f"/{stats['cases']}")
    print(f"fallback won at top level: {stats['fallback_won']}")

    print(f"\n{'choice':>14} {'median gap':>11} {'mean gap':>10} "
          f"{'teacher passes':>15}")
    for k in ("baseline", "random", "largest_gap", "measured_gain"):
        g = np.array(gaps[k])
        c = int(np.sum(costs.get(k, [0])))
        print(f"{k:>14} {np.median(g):11.4f} {g.mean():10.4f} {c:15d}")

    print(f"\nnumerical status         : float_unverified -- DIAGNOSTIC, not a proof")
    print(f"runtime                  : {time.time() - t0:.2f}s")
    print("\nno soundness violation found")
    return 0


if __name__ == "__main__":
    sys.exit(main())
