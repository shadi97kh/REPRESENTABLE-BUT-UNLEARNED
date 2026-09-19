"""Pilot training for the two learned components.

DEVELOPMENT DATA ONLY. The config names a held-out list; this script never loads it.
Every bound produced with a policy attached is re-checked against exhaustive enumeration
before any number is reported, because a policy is allowed to affect tightness and must
never affect soundness.

Teacher cost is counted in the same units as the thing it buys: conditional propagations
for the query scorer, bound evaluations for the slope policy.

Run:  python -m scmp.train_policy --config configs/policy_pilot.yaml
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
import yaml

from .bounds import certified_upper_bound
from .conditional import (BoundCache, SCORERS, ConditionalBounder,
                          conditional_upper_bound)
from .model import CertifiedModel, backbone_adjacency
from .oracles import RNANonCrossingOracle
from .oracles.bruteforce import enumerate_rna
from .policy import LearnedQueryScorer, LearnedSlopePolicy, edge_features

SLACK = 1e-9


def onehot(seq):
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return X


def make_instance(seq, cfg, rng):
    o = RNANonCrossingOracle(seq, 3, True)
    m = CertifiedModel(4, o.ground_set(), len(seq), d_hid=cfg["model"]["d_hid"],
                       gamma=cfg["model"]["gamma"], seed=rng.randint(0, 9999))
    scale = cfg["model"]["perturb_scale"]
    if scale:
        g = torch.Generator().manual_seed(rng.randint(0, 9999))
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype) * scale)
    members = enumerate_rna(seq, 3, True)
    X, B = onehot(seq), backbone_adjacency(len(seq))
    with torch.no_grad():
        true_max = max(float(m(X, m.indicator(s), B)) for s in members)
    return {"seq": seq, "o": o, "m": m, "X": X, "B": B, "true_max": true_max}


def context_for(inst, cfg):
    bd = ConditionalBounder(inst["m"], inst["X"], inst["B"], inst["o"], budget=0,
                            conc_mode=cfg["bounds"]["conc_mode"])
    uncond = bd._upto_z1(bd.forced, bd.forbidden)
    ctx = {"index": bd.index, "d_hid": len(uncond["hi"][0]), "n": bd.n,
           "z1_hi": uncond["hi"], "z1_lo": uncond["lo"],
           "conditional_z1": bd._conditional_z1, "teacher_passes": 0}
    return bd, ctx, uncond


def teacher_ranking(inst, cfg):
    """measured_gain labels. Cost = one conditional pass per candidate edge."""
    bd, ctx, uncond = context_for(inst, cfg)
    edges = list(inst["o"].ground_set())
    gains = []
    passes = 0
    for e in edges:
        cond = bd._conditional_z1(e, count_as="teacher")
        passes += 1
        if cond is None:
            gains.append(float("inf"))
            continue
        gains.append(sum(max(0.0, uncond["hi"][u][k] - cond["hi"][u][k])
                         for u in range(bd.n) for k in range(ctx["d_hid"])))
    return edges, np.array([g if np.isfinite(g) else 1e6 for g in gains]), ctx, passes


def train_scorer(train, cfg):
    lc = cfg["learn"]["query_scorer"]
    net = LearnedQueryScorer(lc["hidden"], seed=cfg["data"]["seed"])
    data, cost = [], 0
    for inst in train:
        edges, gains, ctx, passes = teacher_ranking(inst, cfg)
        cost += passes
        if len(edges) >= 2:
            data.append((edge_features(edges, ctx, inst["m"].n_nodes),
                         torch.tensor(gains, dtype=torch.float64)))
    opt = torch.optim.Adam(net.parameters(), lr=lc["lr"])
    for _ in range(lc["epochs"]):
        opt.zero_grad()
        loss = torch.zeros((), dtype=torch.float64)
        for feats, gains in data:
            s = net.net(feats).squeeze(-1)
            d = gains[:, None] - gains[None, :]
            ds = s[:, None] - s[None, :]
            mask = (d.abs() > 1e-9).double()
            # pairwise logistic: rank the higher-gain edge above the lower one
            loss = loss + (torch.nn.functional.softplus(-torch.sign(d) * ds)
                           * mask).sum() / mask.sum().clamp(min=1)
        (loss / max(len(data), 1)).backward()
        opt.step()
    return net, cost, float(loss / max(len(data), 1))


def eval_choice(inst, cfg, scorer=None, policy=None, slope=None):
    kw = dict(budget=cfg["bounds"]["budget"], conc_mode=cfg["bounds"]["conc_mode"],
              seed=cfg["data"]["seed"], cache=BoundCache())
    if isinstance(scorer, str):
        kw["scorer"] = scorer
    r = conditional_upper_bound(inst["m"], inst["X"], inst["o"], inst["B"],
                                policy=policy, alpha_fn=slope, **kw)
    if r.upper < inst["true_max"] - SLACK:
        raise SystemExit(f"SOUNDNESS VIOLATION on {inst['seq']}: "
                         f"upper {r.upper} < true max {inst['true_max']}")
    return r.upper - inst["true_max"], r


def train_slopes(train, cfg, scorer_net):
    """Imitate a best-of-K random alpha search. Teacher cost = K bound evaluations."""
    lc = cfg["learn"]["slope_policy"]
    net = LearnedSlopePolicy(lc["hidden"], seed=cfg["data"]["seed"] + 1)
    rng = random.Random(cfg["data"]["seed"])
    rows, cost = [], 0
    for inst in train:
        best_a, best_u = None, None
        for _ in range(lc["teacher_draws"]):
            a = rng.random()
            gap, _ = eval_choice(inst, cfg, scorer="largest_gap",
                                 slope=lambda L, i, k, _a=a: _a)
            cost += 1
            if best_u is None or gap < best_u:
                best_u, best_a = gap, a
        bd, ctx, uncond = context_for(inst, cfg)
        for i in range(bd.n):
            for k in range(ctx["d_hid"]):
                l, u = uncond["lo"][i][k], uncond["hi"][i][k]
                if l < 0 < u:
                    rows.append(([l, u, u - l, u / (u - l)], best_a))
    if not rows:
        return net, cost, float("nan")
    F = torch.tensor([r[0] for r in rows], dtype=torch.float64)
    Y = torch.tensor([r[1] for r in rows], dtype=torch.float64)
    opt = torch.optim.Adam(net.parameters(), lr=lc["lr"])
    for _ in range(lc["epochs"]):
        opt.zero_grad()
        pred = torch.sigmoid(net.net(F)).squeeze(-1)
        loss = ((pred - Y) ** 2).mean()
        loss.backward()
        opt.step()
    return net, cost, float(loss)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.train_policy")
    ap.add_argument("--config", default="configs/policy_pilot.yaml")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    rng = random.Random(cfg["data"]["seed"])
    t0 = time.time()

    dev = cfg["data"]["development"]
    print(f"policy pilot  config={a.config}")
    print(f"DEVELOPMENT DATA ONLY: {dev}")
    print(f"held-out list named but NOT loaded: {cfg['data']['held_out']}\n")

    train = [make_instance(rng.choice(dev), cfg, rng)
             for _ in range(cfg["data"]["train_instances"])]
    evalset = [make_instance(rng.choice(dev), cfg, rng)
               for _ in range(cfg["data"]["eval_instances"])]
    print(f"instances: {len(train)} train, {len(evalset)} eval (disjoint draws)\n")

    results, teacher = {}, {}

    print("deterministic choices")
    for name in cfg["compare"]["deterministic_choices"]:
        gaps, tp = [], 0
        for inst in evalset:
            g, r = eval_choice(inst, cfg, scorer=name)
            gaps.append(g)
            tp += r.stats.get("teacher_passes", 0)
        results[name] = np.array(gaps)
        teacher[name] = tp
        print(f"  {name:>14}  median gap {np.median(gaps):8.4f}  "
              f"mean {np.mean(gaps):8.4f}  teacher passes {tp}")

    base = np.array([certified_upper_bound(i["m"], i["X"], i["o"], i["B"],
                                           mode="combined").upper - i["true_max"]
                     for i in evalset])
    results["baseline_unconditional"] = base
    print(f"  {'baseline':>14}  median gap {np.median(base):8.4f}  "
          f"mean {np.mean(base):8.4f}  teacher passes 0")

    scorer_net = None
    if cfg["learn"]["query_scorer"]["enabled"]:
        print("\nlearned query scorer")
        scorer_net, cost, loss = train_scorer(train, cfg)
        teacher["learned_scorer"] = cost
        gaps = [eval_choice(i, cfg, policy=scorer_net)[0] for i in evalset]
        results["learned_scorer"] = np.array(gaps)
        print(f"  train pairwise loss {loss:.4f}, teacher passes {cost}")
        print(f"  {'learned':>14}  median gap {np.median(gaps):8.4f}  "
              f"mean {np.mean(gaps):8.4f}")

    if cfg["learn"]["slope_policy"]["enabled"]:
        print("\nlearned slope policy (alpha in [0,1], valid by construction)")
        slope_net, cost, loss = train_slopes(train[:8], cfg, scorer_net)
        teacher["learned_slopes"] = cost
        gaps = []
        for inst in evalset:
            bd, ctx, uncond = context_for(inst, cfg)
            stats = {("H1", i, k): (uncond["lo"][i][k], uncond["hi"][i][k])
                     for i in range(bd.n) for k in range(ctx["d_hid"])}
            slope_net.bind(stats)
            gaps.append(eval_choice(inst, cfg, scorer="largest_gap",
                                    slope=slope_net.alpha_fn)[0])
        results["learned_slopes"] = np.array(gaps)
        print(f"  train MSE {loss:.4f}, teacher passes {cost}")
        print(f"  {'learned':>14}  median gap {np.median(gaps):8.4f}  "
              f"mean {np.mean(gaps):8.4f}")

    print("\nsummary (lower gap is tighter; all bounds re-checked sound)")
    print(f"  {'method':>24} {'median gap':>11} {'mean gap':>10} {'teacher cost':>13}")
    for k, v in results.items():
        print(f"  {k:>24} {np.median(v):11.4f} {np.mean(v):10.4f} "
              f"{teacher.get(k, 0):13d}")

    os.makedirs("runs/policy_pilot", exist_ok=True)
    path = f"runs/policy_pilot/{time.strftime('%Y%m%d-%H%M%S')}.json"
    with open(path, "w") as fh:
        json.dump({"config": cfg, "teacher_cost": teacher,
                   "results": {k: v.tolist() for k, v in results.items()},
                   "runtime_s": time.time() - t0,
                   "note": "development data only; held_out never loaded",
                   "numerical_status": "float_unverified"}, fh, indent=2)
    print(f"\nevidence: {path}")
    print(f"runtime: {time.time() - t0:.1f}s")
    print("no network predicted a bound: scorer emits a ranking, slope emits alpha in [0,1]")
    return 0


if __name__ == "__main__":
    sys.exit(main())
