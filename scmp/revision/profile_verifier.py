"""P2 profiling: measure where verifier time goes, before optimising anything.

Run this BEFORE any batching or chart-reuse work. The plan is explicit that such
optimisation is justified only where profiling identifies cost, and that outputs
must be verified after each change. This module produces the measurement; it
changes no algorithm.
"""
from __future__ import annotations

import argparse
import json
import os
import time

import numpy as np
import torch

from . import FLOAT_DIAGNOSTIC
from .instrument import COMPONENTS, Meter
from .ladder import LADDER, RULE_HASH, build_case, generate
from .provenance import sha256_file


def _thread_budget(n):
    torch.set_num_threads(int(n))
    for v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS"):
        os.environ[v] = str(int(n))


def profile_case(case, arm, budget, meter=None):
    """Run one arm on one case with every component timed."""
    from ..bounds import certified_upper_bound
    from ..conditional import conditional_upper_bound

    m = meter or Meter()
    with m.timed("chart_compilation"):
        model, X, oracle, B = build_case(case)
        ground = oracle.ground_set()
    with m.timed("support_query", work=1):
        feas = oracle.feasibility()
        fam_size = oracle.count().value

    t_bound = time.perf_counter()
    if arm in ("box", "terminal_only", "intermediate_unconditional"):
        mode = "terminal_only" if arm == "terminal_only" else "combined"
        with m.timed("bound_propagation"):
            r = certified_upper_bound(model, X, oracle, B, mode=mode)
        report = {}
    else:
        with m.timed("bound_propagation"):
            r = conditional_upper_bound(model, X, oracle, B,
                                        budget=budget["conditional_query_budget"],
                                        scorer="largest_gap",
                                        conc_mode="family")
        report = dict(getattr(r, "report", {}) or {})
        if isinstance(getattr(r, "stats", None), dict):
            report.update(r.stats)
    bound_s = time.perf_counter() - t_bound

    with m.timed("incumbent_search", work=1):
        # BoundResult exposes incumbent_value / gap; there is no `lower` field,
        # and getattr(..., "lower", -inf) silently reported every gap as inf.
        incumbent = float(r.incumbent_value)
    m.note_bound(incumbent, float(r.upper),
                 "proved" if float(r.upper) - incumbent <= 1e-9 else "open")

    # Reattribute work reported from inside the bounder. Wrapping
    # conditional_upper_bound in one timer credited 100% to bound_propagation and
    # hid the oracle and cache traffic that is the whole point of conditioning.
    st = dict(getattr(r, "stats", {}) or {})
    ocalls = int(st.get("final_oracle_calls") or st.get("oracle_calls") or 0)
    m.work["support_query"] += ocalls
    m.calls["support_query"] += ocalls
    ch = st.get("cache") or {}
    m.cache["hits"] += int(ch.get("hits", 0))
    m.cache["misses"] += int(ch.get("misses", 0))
    m.calls["cache_access"] += int(ch.get("hits", 0)) + int(ch.get("misses", 0))
    m.work["cache_access"] += int(ch.get("misses", 0))   # only a miss does work
    oracle_stats = oracle.stats() if hasattr(oracle, "stats") else {}
    m.work["chart_compilation"] += int(
        oracle_stats.get("materialized_productions", 0) or 0)

    s = m.summary()
    s["reported_from_bounder"] = {
        "oracle_calls": ocalls,
        "qbounds": st.get("qbounds"),
        "conditional_passes": st.get("conditional_passes"),
        "budget": st.get("budget"),
        "conditional_used": st.get("conditional_used"),
        "fallback_used": st.get("fallback_used"),
        "proved_impossible": st.get("proved_impossible"),
        "relu_unstable": st.get("relu_unstable"),
        "relu_stable_pos": st.get("relu_stable_pos"),
        "relu_stable_neg": st.get("relu_stable_neg"),
        "cache": ch,
        "used": st.get("used"),
    }
    return {"case": {k: case[k] for k in ("rung", "family_tag", "seed", "margin",
                                          "n_candidates", "d_hid")},
            "arm": arm, "family_size": int(fam_size),
            "upper": float(r.upper), "lower": incumbent,
            "gap": float(r.upper) - incumbent,
            "bound_wall_s": round(bound_s, 6),
            "oracle_calls_reported": report.get("oracle_calls"),
            "cache_hits_reported": report.get("cache_hits"),
            "profile": s}


def main(argv=None):
    import yaml
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.profile_verifier",
        description="Component-level profile of the verifier across the frozen "
                    "difficulty ladder. Measures; does not optimise.")
    ap.add_argument("--config", default="configs/revision_profile.yaml")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    cases = generate(cfg["ladder_rungs"], cfg["seeds_per_rung"], cfg["margins"])

    if a.dry_run:
        print(f"dry run\n  ladder rule hash {RULE_HASH}\n  cases {len(cases)}"
              f"\n  arms {cfg['arms']}\n  runs {len(cases) * len(cfg['arms'])}"
              f"\n  threads {cfg['budget']['threads']}, CPU only"
              f"\n  cap {cfg['budget']['max_cpu_core_hours']} core-hours")
        return 0

    _thread_budget(cfg["budget"]["threads"])
    t0 = time.time()
    rows = []
    for case in cases:
        for arm in cfg["arms"]:
            try:
                rows.append(profile_case(case, arm, cfg["budget"]))
            except Exception as exc:
                rows.append({"case": {k: case[k] for k in
                                      ("rung", "family_tag", "seed")},
                             "arm": arm,
                             "error": f"{type(exc).__name__}: {exc}"})

    ok = [r for r in rows if "error" not in r]
    agg = {}
    for arm in cfg["arms"]:
        sub = [r for r in ok if r["arm"] == arm]
        if not sub:
            continue
        tot = {c: float(np.sum([r["profile"]["component_time_s"][c] for r in sub]))
               for c in COMPONENTS}
        s = sum(tot.values()) or 1.0
        agg[arm] = {
            "n": len(sub),
            "component_share": {c: round(tot[c] / s, 4) for c in COMPONENTS},
            "component_time_s": {c: round(tot[c], 4) for c in COMPONENTS},
            "median_bound_wall_s": float(np.median([r["bound_wall_s"] for r in sub])),
            "cache_hit_rate": float(np.mean([r["profile"]["cache_hit_rate"]
                                             for r in sub])),
            "median_gap": float(np.median([r["gap"] for r in sub])),
            "median_oracle_calls": float(np.median(
                [r["profile"]["reported_from_bounder"]["oracle_calls"] for r in sub])),
            "cache_hits_total": int(np.sum(
                [r["profile"]["cache"]["hits"] for r in sub])),
            "cache_misses_total": int(np.sum(
                [r["profile"]["cache"]["misses"] for r in sub])),
            "fallback_used_total": int(np.sum(
                [(r["profile"]["reported_from_bounder"]["fallback_used"] or 0)
                 for r in sub])),
            "median_unstable_relus": float(np.median(
                [(r["profile"]["reported_from_bounder"]["relu_unstable"] or 0)
                 for r in sub])),
        }

    print(f"ladder rule hash {RULE_HASH}   cases {len(cases)}   "
          f"runs {len(rows)}   errors {len(rows) - len(ok)}\n")
    print(f"{'arm':<28}{'n':>4}{'med s':>9}{'med gap':>11}{'oracle':>9}"
          f"{'cache h/m':>14}{'unstab':>8}")
    for arm, v in agg.items():
        print(f"{arm:<28}{v['n']:>4}{v['median_bound_wall_s']:>9.4f}"
              f"{v['median_gap']:>11.2f}{v['median_oracle_calls']:>9.0f}"
              f"{v['cache_hits_total']:>7}/{v['cache_misses_total']:<6}"
              f"{v['median_unstable_relus']:>8.0f}")
    print("\ncomponent time share")
    for arm, v in agg.items():
        top = sorted(v["component_share"].items(), key=lambda kv: -kv[1])[:3]
        print(f"  {arm:<28}" + "  ".join(f"{k}={p:.0%}" for k, p in top if p > 0))

    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir,
                        f"profile_verifier-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump({"kind": "revision_profile_verifier", "config": a.config,
                   "config_sha256": sha256_file(a.config),
                   "ladder_rule_hash": RULE_HASH, "ladder": LADDER,
                   "evidence_kind": FLOAT_DIAGNOSTIC,
                   "aggregate": agg, "rows": rows,
                   "cpu_core_hours": round((time.time() - t0)
                                           * cfg["budget"]["threads"] / 3600, 5),
                   "runtime_s": round(time.time() - t0, 2)}, fh, indent=1,
                  default=str)
    with open(path + ".sha256", "w") as fh:
        fh.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
