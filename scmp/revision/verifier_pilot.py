"""P2: verifier arms under protocol v2, on the frozen difficulty ladder.

WHY PROTOCOL V1'S METRIC IS REPLACED. v1's primary quantity was a median relative
gap reduction. On the v1 families the plain box baseline already closed 99.1% of
instances to gap exactly zero, so the metric was 0/0 and G1/G2 were recorded
UNDETERMINED. Those gates are retained historically and are NOT re-scored here.
v2 declares three quantities instead:

  decision coverage        fraction of instances where an arm CERTIFIES the
                           declared decision within the shared budget. A decision
                           is settled either way -- proving the optimum is below
                           tau, or exhibiting a feasible point at or above it.
  capped time-to-resolution first time the gap reaches tolerance; unresolved
                           instances are recorded AT the cap and counted.
  anytime gap AUC          (1/T) integral of min(1, gap(t)/scale) dt on a scale
                           shared by every arm on that instance.

WHAT IS NOT DONE HERE. No case is added, dropped or reweighted after an arm has
run on it; the generation rule is hashed in `ladder.py` before any case is drawn.
Rung 0 is the preserved v1 easy set and always appears in the aggregate.
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import time

import numpy as np
import torch

from . import FLOAT_DIAGNOSTIC
from .instrument import Meter
from .ladder import LADDER, RULE_HASH, build_case, generate
from .provenance import sha256_file

ARMS = ("box", "terminal_only", "intermediate_unconditional", "fixed_conditional",
        "deterministic_selective", "random_selective", "lp_polytope_external")


# ------------------------------------------------------- external references
def exhaustive_optimum(model, X, oracle, B, max_size, max_edges=80,
                       node_budget=200000):
    """Exact max over the family by enumeration. Ground truth, not a bound.

    OPTIONAL reference. Where it is too expensive it returns None with a reason,
    and the case is reported WITHOUT ground truth rather than with a guess.

    Cost control matters here. `feasibility()` runs the family's full O(n^3) DP,
    and calling it at every DFS node made a single 12-nt RNA case exceed two
    minutes. Three guards: memoise feasibility on the prefix, cap the candidate
    edges (the DFS is exponential in them, not in the family size), and cap the
    number of DFS nodes so a bad case degrades to "unavailable" instead of
    hanging the run.
    """
    edges = list(oracle.ground_set())
    if len(edges) > max_edges:
        return None, {"reason": f"{len(edges)} candidate edges exceeds cap "
                                f"{max_edges}; ground truth not attempted"}
    n = oracle.count().value
    if n > max_size:
        return None, {"reason": f"family size {n} exceeds cap {max_size}"}

    feas_cache = {}

    def feasible(prefix):
        key = prefix
        v = feas_cache.get(key)
        if v is None:
            v = oracle.feasibility(forced=key).feasible
            feas_cache[key] = v
        return v

    eset = set(edges)

    def is_member(prefix):
        """Exactly one structure satisfies 'these edges present, all others absent'."""
        rest = tuple(sorted(eset - set(prefix)))
        return oracle.count(forced=tuple(sorted(prefix)),
                            forbidden=rest).value == 1

    best, best_z, seen, nodes = float("-inf"), None, 0, 0

    def dfs(k, prefix):
        nonlocal best, best_z, seen, nodes
        nodes += 1
        if nodes > node_budget:
            raise RuntimeError(f"DFS exceeded {node_budget} nodes")
        if k == len(edges):
            # `feasible(S)` means "extendable to a member", which is NOT the same
            # as "is a member" unless the family is downward-closed. Layered-DAG
            # paths are not: a partial path is extendable but is not itself a
            # path. Evaluating those inflated the maximum and made every sound
            # bound look unsound. Membership is tested by pinning the complement.
            if not is_member(prefix):
                return
            seen += 1
            with torch.no_grad():
                v = float(model(X, model.indicator(prefix), B))
            if v > best:
                best, best_z = v, prefix
            return
        dfs(k + 1, prefix)                                   # exclude edge k
        ext = prefix + (edges[k],)
        if feasible(ext):                                    # prune infeasible
            dfs(k + 1, ext)

    try:
        dfs(0, ())
    except RuntimeError as exc:
        return None, {"reason": str(exc)}
    if best_z is None:
        return None, {"reason": "no feasible member found"}
    return best, {"members_evaluated": seen, "dfs_nodes": nodes,
                  "feasibility_cache_entries": len(feas_cache),
                  "argmax": [list(e) for e in best_z],
                  "family_count_from_oracle": int(n)}


def lp_polytope_bound(model, X, oracle, B):
    """Independent bound: LP over an outer relaxation of the family polytope.

    max a^T z  s.t.  0 <= z <= 1 and, for every node, the incident z sum to <= 1.
    Degree-at-most-one is valid for all three families here (matchings, layered
    paths, non-crossing pairings), so this is a genuine outer relaxation -- looser
    than the exact oracle by construction, which is the point: it is a DIFFERENT
    bounding strategy, not ours.

    The residual is bounded by its box envelope, exactly as `terminal_only` does,
    because the LP does not model the network.
    """
    from scipy.optimize import linprog

    from ..bounds import ResidualBounder
    edges = list(oracle.ground_set())
    if not edges:
        return None, {"reason": "empty ground set"}
    a = model.additive_coefficients(X)
    c = -np.array([a[e] for e in edges])              # linprog minimises
    nodes = sorted({v for e in edges for v in e})
    A_ub = np.zeros((len(nodes), len(edges)))
    for r, v in enumerate(nodes):
        for k, e in enumerate(edges):
            if v in e:
                A_ub[r, k] = 1.0
    res = linprog(c, A_ub=A_ub, b_ub=np.ones(len(nodes)),
                  bounds=[(0, 1)] * len(edges), method="highs")
    if not res.success:
        return None, {"reason": f"linprog failed: {res.message}"}
    lin = model.anchor_bias(X) + float(-res.fun)
    extra = 0.0
    if model.residual_active:
        lo, hi, _ = ResidualBounder(model, X, B).propagate()
        use = hi if model.gamma >= 0 else lo
        extra = model.gamma * (use.const + float(np.sum(np.maximum(use.coef, 0.0))))
    return lin + extra, {"lp_status": res.status, "n_constraints": len(nodes),
                         "relaxation": "degree<=1 outer relaxation; ignores "
                                       "non-crossing and layer structure"}


# ---------------------------------------------------- deterministic selective
def selective_decision(model, X, oracle, B, cfg_sel, meter):
    """Decide whether to condition at all, and record the action taken.

    Three actions, and 'skip' is a first-class one. A policy that is REQUIRED to
    condition wastes oracle calls on instances that cheap refinement already
    resolves -- which is what the v1 conditional arms did, and why they lost at
    matched wall time despite producing tighter bounds per query.
    """
    from ..bounds import certified_upper_bound
    with meter.timed("bound_propagation"):
        cheap = certified_upper_bound(model, X, oracle, B, mode="combined")
    gap = float(cheap.gap)
    meter.note_bound(float(cheap.incumbent_value), float(cheap.upper))

    st = dict(cheap.stats or {})
    unst = int(st.get("relu_unstable", 0) or 0)
    stable = int(st.get("relu_stable_pos", 0) or 0) + int(st.get("relu_stable_neg", 0) or 0)
    frac_stable = stable / max(stable + unst, 1)

    if gap <= cfg_sel["skip_if_gap_below"]:
        return "skip_already_resolved", cheap, {"gap": gap, "frac_stable": frac_stable}
    if frac_stable >= cfg_sel["skip_if_stable_fraction"]:
        return "skip_few_unstable", cheap, {"gap": gap, "frac_stable": frac_stable}
    return "condition", cheap, {"gap": gap, "frac_stable": frac_stable,
                                "unstable": unst}


# ------------------------------------------------------------------ one arm
def run_arm(arm, case, cfg, scale_hint=None):
    from ..bounds import certified_upper_bound
    from ..conditional import conditional_upper_bound

    m = Meter()
    with m.timed("chart_compilation"):
        model, X, oracle, B = build_case(case)
    action, extra = None, {}

    if arm == "lp_polytope_external":
        with m.timed("bound_propagation"):
            ub, info = lp_polytope_bound(model, X, oracle, B)
        with m.timed("incumbent_search", work=1):
            base = certified_upper_bound(model, X, oracle, B, mode="combined")
            inc = float(base.incumbent_value)
        if ub is None:
            return {"arm": arm, "error": info.get("reason"), "case": case}
        upper, extra = float(ub), info
        m.note_bound(inc, upper)
        st = {}
    elif arm in ("box", "intermediate_unconditional", "terminal_only"):
        mode = "terminal_only" if arm == "terminal_only" else "combined"
        with m.timed("bound_propagation"):
            r = certified_upper_bound(model, X, oracle, B, mode=mode)
        upper, inc, st = float(r.upper), float(r.incumbent_value), dict(r.stats or {})
        m.note_bound(inc, upper)
    elif arm in ("fixed_conditional", "random_selective"):
        scorer = "random" if arm == "random_selective" else "largest_gap"
        with m.timed("bound_propagation"):
            r = conditional_upper_bound(model, X, oracle, B,
                                        budget=cfg["selective"]["budget"],
                                        scorer=scorer, seed=case["seed"],
                                        conc_mode="family")
        upper, inc, st = float(r.upper), float(r.incumbent_value), dict(r.stats or {})
        m.note_bound(inc, upper)
        action = "condition"
    elif arm == "deterministic_selective":
        action, cheap, extra = selective_decision(model, X, oracle, B,
                                                  cfg["selective"], m)
        if action == "condition":
            with m.timed("bound_propagation"):
                r = conditional_upper_bound(model, X, oracle, B,
                                            budget=cfg["selective"]["budget"],
                                            scorer="largest_gap", seed=case["seed"],
                                            conc_mode="family")
            # retain the MINIMUM of valid bounds: the cheap bound stays valid
            if float(r.upper) <= float(cheap.upper):
                upper, inc, st = float(r.upper), float(r.incumbent_value), dict(r.stats or {})
            else:
                upper, inc, st = float(cheap.upper), float(cheap.incumbent_value), dict(cheap.stats or {})
                extra["fallback_to_cheap"] = True
            m.note_bound(inc, upper)
        else:
            upper, inc, st = float(cheap.upper), float(cheap.incumbent_value), dict(cheap.stats or {})
    else:
        raise ValueError(arm)

    # oracle / cache work reported from inside the bounder
    ocalls = int(st.get("final_oracle_calls") or st.get("oracle_calls") or 0)
    m.work["support_query"] += ocalls
    ch = st.get("cache") or {}
    m.cache["hits"] += int(ch.get("hits", 0))
    m.cache["misses"] += int(ch.get("misses", 0))

    gap = upper - inc
    return {"arm": arm, "case": {k: case[k] for k in
                                 ("rung", "rung_name", "family_tag", "seed",
                                  "margin", "n_candidates", "d_hid")},
            "upper": upper, "incumbent": inc, "gap": gap,
            "negative_gap_defect": bool(gap < 0),
            "action": action, "action_detail": extra,
            "oracle_calls": ocalls,
            "cache": {"hits": int(ch.get("hits", 0)),
                      "misses": int(ch.get("misses", 0))},
            "fallback_used": st.get("fallback_used"),
            "relu_unstable": st.get("relu_unstable"),
            "profile": m.summary(), "meter": m}


# ------------------------------------------------------- protocol v2 metrics
def score_case(rows_by_arm, exact, cap_s, floor):
    """Apply the three declared v2 quantities to one case, across all arms.

    `tau` is derived from the case margin and the SHARED incumbent, so every arm
    faces the same decision. The AUC scale is likewise shared, taken from the box
    arm's initial gap, so the AUCs are comparable rather than each arm being
    graded on its own scale.
    """
    ref = rows_by_arm.get("box") or next(iter(rows_by_arm.values()))
    inc = ref["incumbent"]
    margin = ref["case"]["margin"]
    tau = inc + margin * abs(inc)
    scale = max(ref["gap"], floor)

    out = {}
    for arm, r in rows_by_arm.items():
        if "error" in r:
            out[arm] = {"error": r["error"]}
            continue
        m = r["meter"]
        # decision: is max f over the family below tau?
        certified_below = r["upper"] < tau
        certified_at_or_above = r["incumbent"] >= tau
        auc, negatives = m.anytime_auc(cap_s, scale)
        ttr = m.first_resolution_s
        resolved = ttr is not None
        out[arm] = {
            "upper": r["upper"], "incumbent": r["incumbent"], "gap": r["gap"],
            "tau": tau, "scale": scale,
            "decision_certified": bool(certified_below or certified_at_or_above),
            "decision_direction": ("below" if certified_below else
                                   "at_or_above" if certified_at_or_above else None),
            "capped_time_to_resolution_s": ttr if resolved else cap_s,
            "resolved": bool(resolved),
            "unresolved_reason": None if resolved else "gap above tolerance at cap",
            "anytime_gap_auc": auc,
            "negative_gap_points": negatives,
            "negative_gap_defect": r["negative_gap_defect"],
            "wall_s": m.summary()["wall_s"],
            "oracle_calls": r["oracle_calls"],
            "action": r["action"],
            "sound_vs_exact": (None if exact is None else bool(r["upper"] >= exact - 1e-6)),
            "true_gap_vs_exact": (None if exact is None else r["upper"] - exact),
        }
    return out, {"tau": tau, "scale": scale, "shared_incumbent": inc,
                 "exact_optimum": exact}


def main(argv=None):
    import yaml
    ap = argparse.ArgumentParser(
        prog="python -m scmp.revision.verifier_pilot",
        description="Verifier arms on the frozen difficulty ladder under "
                    "protocol v2. v1 gates are retained, never re-scored.")
    ap.add_argument("--config", default="configs/revision_verifier.yaml")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--out-dir", default="runs/revision")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    cases = generate(cfg["ladder_rungs"], cfg["seeds_per_rung"], cfg["margins"])
    arms = cfg["arms"]

    if a.dry_run:
        print("dry run -- nothing verified\n")
        print(f"ladder rule hash   {RULE_HASH}  (frozen before any case drawn)")
        print(f"cases              {len(cases)} over rungs {cfg['ladder_rungs']}")
        print(f"arms               {len(arms)}: {', '.join(arms)}")
        print(f"runs               {len(cases) * len(arms)}")
        print(f"budget             {cfg['budget']['time_limit_s']}s/instance, "
              f"{cfg['budget']['threads']} thread(s), CPU only, "
              f"cap {cfg['budget']['max_cpu_core_hours']} core-hours")
        print(f"primary metric     {cfg['protocol_v2']['primary']}")
        print(f"external reference exhaustive enumeration (<= "
              f"{cfg['external_reference']['exhaustive_enumeration']['max_family_size']}"
              f" members) + LP polytope relaxation")
        print("unavailable        " + ", ".join(
            u["name"] for u in cfg["external_reference"]["unavailable"]))
        print("\nrung 0 is the preserved v1 easy set and stays in the aggregate.")
        print("v1 gates are retained unchanged and are not re-scored.")
        return 0

    torch.set_num_threads(int(cfg["budget"]["threads"]))
    t0 = time.time()
    cap = cfg["budget"]["time_limit_s"]
    floor = float(cfg["protocol_v2"]["secondary"][1]["anytime_gap_auc"]
                  ["task_scale_floor"])
    maxfam = cfg["external_reference"]["exhaustive_enumeration"]["max_family_size"]

    scored, raw = [], []
    for ci, case in enumerate(cases):
        by_arm = {}
        for arm in arms:
            try:
                by_arm[arm] = run_arm(arm, case, cfg)
            except Exception as exc:
                by_arm[arm] = {"arm": arm, "error": f"{type(exc).__name__}: {exc}",
                               "case": case}
        exact, exinfo = (None, {})
        if cfg["external_reference"]["exhaustive_enumeration"]["use"]:
            try:
                model, X, oracle, B = build_case(case)
                exact, exinfo = exhaustive_optimum(model, X, oracle, B, maxfam)
            except Exception as exc:
                exinfo = {"reason": f"{type(exc).__name__}: {exc}"}
        s, meta = score_case(by_arm, exact, cap, floor)
        scored.append({"case": case, "meta": meta, "exact_info": exinfo, "arms": s})
        raw.append({k: {kk: vv for kk, vv in v.items() if kk != "meter"}
                    for k, v in by_arm.items()})

    # ---------------- aggregate, per rung AND overall (rung 0 never dropped) --
    def agg(sel, label):
        out = {}
        for arm in arms:
            vals = [c["arms"][arm] for c in sel
                    if arm in c["arms"] and "error" not in c["arms"][arm]]
            if not vals:
                continue
            out[arm] = {
                "n": len(vals),
                "decision_coverage": float(np.mean([v["decision_certified"] for v in vals])),
                "median_time_to_resolution_s": float(np.median(
                    [v["capped_time_to_resolution_s"] for v in vals])),
                "resolved_fraction": float(np.mean([v["resolved"] for v in vals])),
                "mean_anytime_auc": float(np.mean([v["anytime_gap_auc"] for v in vals])),
                "median_gap": float(np.median([v["gap"] for v in vals])),
                "median_wall_s": float(np.median([v["wall_s"] for v in vals])),
                "median_oracle_calls": float(np.median([v["oracle_calls"] for v in vals])),
                "unresolved": int(sum(1 for v in vals if not v["resolved"])),
                "negative_gap_defects": int(sum(1 for v in vals if v["negative_gap_defect"])),
                "unsound_vs_exact": int(sum(1 for v in vals
                                            if v["sound_vs_exact"] is False)),
                "median_true_gap": (float(np.median(
                    [v["true_gap_vs_exact"] for v in vals
                     if v["true_gap_vs_exact"] is not None]))
                    if any(v["true_gap_vs_exact"] is not None for v in vals) else None),
            }
        return {"label": label, "arms": out}

    aggregates = {"overall": agg(scored, "overall (rung 0 included)")}
    for r in sorted({c["case"]["rung"] for c in scored}):
        aggregates[f"rung{r}"] = agg([c for c in scored if c["case"]["rung"] == r],
                                     f"rung {r}")

    # ---------------- selective action census -------------------------------
    acts = {}
    for c in scored:
        v = c["arms"].get("deterministic_selective", {})
        acts[v.get("action")] = acts.get(v.get("action"), 0) + 1

    # ---------------- policy gate -------------------------------------------
    ov = aggregates["overall"]["arms"]
    det = ov.get("deterministic_selective", {})
    best_cheap = min((ov[a] for a in ("box", "intermediate_unconditional")
                      if a in ov), key=lambda d: d["mean_anytime_auc"], default=None)
    gate = {"criterion": "deterministic selective must beat the best cheap arm on "
                         "decision coverage OR anytime AUC at equal wall budget",
            "deterministic": {k: det.get(k) for k in
                              ("decision_coverage", "mean_anytime_auc",
                               "median_wall_s")},
            "best_cheap": ({k: best_cheap.get(k) for k in
                            ("decision_coverage", "mean_anytime_auc",
                             "median_wall_s")} if best_cheap else None)}
    if det and best_cheap:
        gate["useful_region"] = bool(
            det["decision_coverage"] > best_cheap["decision_coverage"]
            or det["mean_anytime_auc"] < best_cheap["mean_anytime_auc"])
    else:
        gate["useful_region"] = False
    gate["policy_trained"] = False
    gate["reason"] = ("policy NOT trained: the plan gates it on deterministic "
                      "selective conditioning showing a useful cost-benefit "
                      "region, and it does not"
                      if not gate["useful_region"] else
                      "deterministic selective shows a useful region; a policy "
                      "may now be trained under P2's stated constraints")

    # ---------------- report -------------------------------------------------
    print(f"ladder rule hash {RULE_HASH}   cases {len(cases)}   "
          f"runs {len(cases) * len(arms)}\n")
    for key in ["overall"] + [f"rung{r}" for r in sorted(
            {c['case']['rung'] for c in scored})]:
        A = aggregates[key]["arms"]
        print(f"{aggregates[key]['label']}")
        print(f"  {'arm':<28}{'cov':>7}{'AUC':>8}{'ttr s':>9}{'unres':>7}"
              f"{'wall s':>9}{'oracle':>8}{'true gap':>10}")
        for arm in arms:
            v = A.get(arm)
            if not v:
                continue
            tg = "-" if v["median_true_gap"] is None else f"{v['median_true_gap']:.3f}"
            print(f"  {arm:<28}{v['decision_coverage']:>7.1%}"
                  f"{v['mean_anytime_auc']:>8.3f}"
                  f"{v['median_time_to_resolution_s']:>9.3f}{v['unresolved']:>7}"
                  f"{v['median_wall_s']:>9.4f}{v['median_oracle_calls']:>8.0f}"
                  f"{tg:>10}")
        print()

    print(f"selective actions: {acts}")
    unsound = sum(v["unsound_vs_exact"] for A in aggregates.values()
                  for v in A["arms"].values())
    negs = sum(v["negative_gap_defects"] for A in aggregates.values()
               for v in A["arms"].values())
    print(f"soundness vs exact optimum: {unsound} violation(s); "
          f"negative-gap defects: {negs}")
    print(f"\npolicy gate: useful_region={gate['useful_region']} -> "
          f"{gate['reason']}")

    res = {"kind": "revision_verifier_pilot", "config": a.config,
           "config_sha256": sha256_file(a.config), "ladder_rule_hash": RULE_HASH,
           "ladder": LADDER, "protocol_v2": cfg["protocol_v2"],
           "evidence_kind": FLOAT_DIAGNOSTIC,
           "v1_gates": "retained unchanged; not re-scored under v2",
           "external_reference": cfg["external_reference"],
           "aggregates": aggregates, "selective_actions": acts,
           "policy_gate": gate, "cases": scored, "raw": raw,
           "runtime_s": round(time.time() - t0, 2),
           "cpu_core_hours": round((time.time() - t0)
                                   * cfg["budget"]["threads"] / 3600, 5)}
    os.makedirs(a.out_dir, exist_ok=True)
    path = os.path.join(a.out_dir,
                        f"verifier_pilot-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump(res, fh, indent=1, default=str)
    with open(path + ".sha256", "w") as fh:
        fh.write(sha256_file(path) + "  " + os.path.basename(path) + "\n")
    print(f"\nwritten {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
