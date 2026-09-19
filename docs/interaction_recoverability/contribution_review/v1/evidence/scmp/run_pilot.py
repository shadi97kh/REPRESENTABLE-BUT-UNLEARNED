"""Capped pilot runner.

REFUSES TO RUN without `--dry-run` unless `pilot.authorized_budget` is set.

Live mode logs one record per (family, instance, seed, arm) with wall time, oracle
calls, timeout status and the hardware setting in force, so nothing in the results table
can be traced back to an unlogged choice.

Verifier arms are compared at MATCHED WALL TIME: every arm gets the same per-instance
time limit inside the refinement loop, so a tighter bound bought with more time is not
counted as a tighter method.

Run:  python -m scmp.run_pilot --config configs/pilot.yaml [--dry-run]
"""
from __future__ import annotations

import argparse
import json
import os
import platform
import random
import sys
import time

import numpy as np
import torch
import yaml

from .bounds import certified_upper_bound
from .conditional import BoundCache, conditional_upper_bound
from .model import CertifiedModel, backbone_adjacency
from .oracles import (EdgeRelabeledOracle, LayeredDAGPathOracle, PathMatchingOracle,
                      RNANonCrossingOracle)
from .refine import refine
from .validate_protocol import validate

ARMS = ("box", "terminal_only", "intermediate_unconditional",
        "fixed_conditional", "learned_conditional")


# ------------------------------------------------------------------ hardware
def reserve_hardware(pilot):
    hw = pilot.get("hardware", {}) or {}
    info = {"platform": platform.platform(), "python": platform.python_version(),
            "torch": torch.__version__, "device": "cpu", "gpu_index": None,
            "gpu_memory_fraction": None, "gpu_name": None,
            "requested_gpu": bool(hw.get("use_gpu"))}
    if not hw.get("use_gpu"):
        return info
    from .gpu import reserve
    idx = reserve(fraction=hw.get("gpu_memory_fraction", 0.15), verbose=False)
    if idx is None:
        info["note"] = "no device with enough free memory; running on CPU"
        return info
    info.update(device="cuda:0", gpu_index=idx,
                gpu_memory_fraction=hw.get("gpu_memory_fraction", 0.15),
                gpu_name=torch.cuda.get_device_name(0))
    return info


# ------------------------------------------------------------------ families
def build_family(spec, rng):
    """Returns (oracle_with_nodepair_ground_set, n_nodes)."""
    name = spec["name"]
    if name == "layered_dag_path":
        sizes = [1] + [rng.randint(1, spec["width"]) for _ in range(spec["n_layers"] - 1)]
        edges = [(l, u, v) for l in range(len(sizes) - 1)
                 for u in range(sizes[l]) for v in range(sizes[l + 1])
                 if rng.random() < 0.75]
        base = LayeredDAGPathOracle(sizes, edges)
        flat, i = {}, 0
        for li, s in enumerate(sizes):
            for u in range(s):
                flat[(li, u)] = i; i += 1
        mapping = {e: tuple(sorted((flat[(e[0], e[1])], flat[(e[0] + 1, e[2])])))
                   for e in base.ground_set()}
        if not mapping:
            return None, 0
        return EdgeRelabeledOracle(base, mapping), i
    if name == "path_matching":
        n = spec["n_nodes"]
        return PathMatchingOracle(n), n
    if name == "rna_noncrossing":
        seq = rng.choice(spec["sequences"])
        return RNANonCrossingOracle(seq, 3, True), len(seq)
    raise ValueError(f"unknown family {name!r}")


def arm_bound_fn(arm, m, X, o, B, budget):
    if arm == "box":
        return lambda f, b: certified_upper_bound(m, X, o, B, f, b, mode="combined")
    if arm == "terminal_only":
        return lambda f, b: certified_upper_bound(m, X, o, B, f, b, mode="terminal_only")
    if arm == "intermediate_unconditional":
        return lambda f, b: certified_upper_bound(m, X, o, B, f, b, mode="combined")
    if arm == "fixed_conditional":
        return lambda f, b: conditional_upper_bound(m, X, o, B, f, b, budget=budget,
                                                    scorer="largest_gap",
                                                    cache=BoundCache())
    if arm == "learned_conditional":
        return lambda f, b: conditional_upper_bound(m, X, o, B, f, b, budget=budget,
                                                    scorer="measured_gain",
                                                    cache=BoundCache())
    raise ValueError(arm)


def run_verifier(cfg, hw, rng, log):
    matched = float(cfg.get("matched_wall_time_s", 1.0))
    budget = cfg["budget"]
    for spec in cfg["families"]:
        for inst in range(spec["instances"]):
            o, n_nodes = build_family(spec, rng)
            if o is None or not o.ground_set():
                log.append({"kind": "void", "family": spec["name"], "instance": inst,
                            "reason": "empty ground set"})
                continue
            for seed in cfg["seeds"]:
                m = CertifiedModel(4, o.ground_set(), n_nodes, d_hid=4, gamma=1.0,
                                   seed=seed)
                g = torch.Generator().manual_seed(seed + 1000 * inst)
                with torch.no_grad():
                    for p in m.parameters():
                        p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype))
                X = torch.rand(n_nodes, 4, dtype=torch.float64,
                               generator=torch.Generator().manual_seed(seed))
                B = backbone_adjacency(n_nodes)
                cnt = o.count().value
                if cnt <= 1:
                    log.append({"kind": "void", "family": spec["name"],
                                "instance": inst, "seed": seed,
                                "reason": "empty or singleton family"})
                    continue
                for arm in cfg["verifier_arms"]:
                    t0 = time.time()
                    try:
                        r = refine(m, X, o, B, max_nodes=cfg["max_nodes"],
                                   time_limit=matched,
                                   bound_fn=arm_bound_fn(arm, m, X, o, B, budget))
                        wall = time.time() - t0
                        denom = max(abs(r.incumbent), 1.0)
                        log.append({"kind": "verifier", "family": spec["name"],
                                    "instance": inst, "seed": seed, "arm": arm,
                                    "gap_rel": (r.global_upper - r.incumbent) / denom,
                                    "upper": r.global_upper,
                                    "incumbent": r.incumbent,
                                    "status": r.status,
                                    "timeout": r.status == "unresolved",
                                    "wall_s": wall, "matched_wall_s": matched,
                                    "nodes": r.stats.get("expanded", 0),
                                    "bounds": r.stats.get("bounds", 0),
                                    "budget": budget, "family_size": cnt,
                                    "numerical_status": r.numerical_status,
                                    "device": hw["device"]})
                    except Exception as exc:
                        log.append({"kind": "error", "family": spec["name"],
                                    "instance": inst, "seed": seed, "arm": arm,
                                    "error": f"{type(exc).__name__}: {exc}"})


# ------------------------------------------------------------------- predictor
def run_predictor(cfg, hw, log):
    """Gene-disjoint RNA dev split. Marginal-label fit only (preregistration section 6)."""
    pc = cfg.get("predictor", {})
    if not pc.get("enabled"):
        return
    from .oracles.bruteforce import enumerate_rna
    sys.path.insert(0, ".")
    from experiments._huesken import build_graphs, gene_map, load_verified
    from certmp.ensemble import bpp_lattice
    from certmp import data as cdata

    records, rep, ver = load_verified()
    gm = gene_map()
    rng = random.Random(cfg["seeds"][0])
    subset = records[:pc["max_sequences"]]
    seqs = [r["sequence"] for r in subset]
    y = np.array([r["efficacy"] for r in subset], dtype=float)
    genes = np.array([gm.get(s, "?") for s in seqs])
    uniq = sorted(set(genes) - {"?"})
    rng.shuffle(uniq)
    held = set(uniq[:max(1, len(uniq) // 4)])
    te = np.where(np.isin(genes, list(held)))[0]
    tr = np.where(~np.isin(genes, list(held)))[0]
    log.append({"kind": "predictor_split", "n": len(subset), "train": len(tr),
                "test": len(te), "held_out_genes": sorted(held),
                "unit": "gene", "note": "development data; marginal label"})
    if len(te) < 10 or len(tr) < 20:
        log.append({"kind": "predictor_abstain", "reason": "split too small"})
        return

    device = torch.device(hw["device"])

    def positional(seq):
        idx = {c: t for t, c in enumerate("ACGU")}
        L = len(seq)
        Xs = np.zeros((L, L * 4))
        for t, c in enumerate(seq.replace("T", "U")):
            if c in idx:
                Xs[t, t * 4 + idx[c]] = 1.0
        return Xs

    graphs = build_graphs(seqs, positional)
    Xall = torch.tensor(np.stack([g["X"] for g in graphs]), dtype=torch.float64)
    Aall = torch.tensor(np.stack([g["A_mfe"] for g in graphs]), dtype=torch.float64)
    Yall = torch.tensor(y, dtype=torch.float64)

    device = torch.device(hw["device"])

    class PlainMPNN(torch.nn.Module):
        """Unconstrained baseline: no anchor, no residual decomposition, no certificate.
        Same depth, width and readout as the certifiable model, so the comparison
        isolates the decomposition rather than capacity."""
        def __init__(self, d_in, d_hid):
            super().__init__()
            self.l0 = torch.nn.Linear(d_in, d_hid, dtype=torch.float64)
            self.l1 = torch.nn.Linear(d_hid, d_hid, dtype=torch.float64)
            self.ro = torch.nn.Linear(d_hid, 1, dtype=torch.float64)

        def forward(self, Xb, Ab):
            H = torch.relu(torch.bmm(Ab, self.l0(Xb)))
            H = torch.relu(torch.bmm(Ab, self.l1(H)))
            return self.ro(H.sum(1)).squeeze(-1)

    n_nodes = Xall.shape[1]
    d_in = Xall.shape[2]
    tri = torch.tensor(tr, device=device)
    tei = torch.tensor(te, device=device)
    Xd, Ad, Yd = Xall.to(device), Aall.to(device), Yall.to(device)

    for arm in pc["arms"]:
        for seed in cfg["seeds"]:
            t0 = time.time()
            torch.manual_seed(seed)
            if arm == "unconstrained":
                net = PlainMPNN(d_in, pc["d_hid"]).to(device)
                params = list(net.parameters())
                fwd = lambda idx: net(Xd[idx], Ad[idx])
                certifiable, gamma = False, None
            else:
                gamma = 1.0 if arm.endswith("gamma1") else 0.0
                cm = CertifiedModel(d_in, tuple(), n_nodes, d_hid=pc["d_hid"],
                                    gamma=gamma, seed=seed).to(device)
                zero_z = torch.zeros(0, dtype=torch.float64, device=device)
                params = list(cm.parameters())

                def fwd(idx, cm=cm, zero_z=zero_z):
                    return torch.stack([cm(Xd[i], zero_z, Ad[i]) for i in idx])
                certifiable = True
            opt = torch.optim.Adam(params, lr=pc["lr"])
            for _ in range(pc["epochs"]):
                opt.zero_grad()
                loss = ((fwd(tri) - Yd[tri]) ** 2).mean()
                loss.backward()
                opt.step()
            with torch.no_grad():
                pv = fwd(tei).detach().cpu().numpy()
            log.append({"kind": "predictor", "arm": arm, "seed": seed,
                        "certifiable": certifiable, "gamma": gamma,
                        "spearman": cdata.spearman(pv, y[te]),
                        "pearson": cdata.pearson(pv, y[te]),
                        "train_loss": float(loss.detach()), "wall_s": time.time() - t0,
                        "device": hw["device"], "epochs": pc["epochs"],
                        "n_train": len(tr), "n_test": len(te),
                        "label_semantics": "marginal over latent structures"})


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.run_pilot")
    ap.add_argument("--config", default="configs/pilot.yaml")
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args(argv)
    cfg = yaml.safe_load(open(a.config))
    pilot = cfg["pilot"]
    authorized = pilot.get("authorized_budget")

    proto = yaml.safe_load(open(pilot["protocol"]))
    rep = validate(proto)
    print(f"pilot          : {a.config}")
    print(f"protocol       : {pilot['protocol']}  "
          f"({'accepted' if not rep.errors else 'REJECTED'}, {rep.checks} obligations)")
    if rep.errors:
        for e in rep.errors:
            print(f"  ERROR {e}")
        return 2
    if not a.dry_run and authorized is None:
        print("\nREFUSING TO RUN: no authorised compute budget "
              "(`pilot.authorized_budget` is null).")
        return 3

    hw = reserve_hardware(pilot) if not a.dry_run else {"device": "cpu",
                                                        "gpu_index": None}
    print(f"authorised     : {authorized}")
    print(f"hardware       : device={hw['device']} gpu_index={hw.get('gpu_index')} "
          f"{hw.get('gpu_name') or ''} frac={hw.get('gpu_memory_fraction')}")
    print(f"mode           : {'DRY RUN' if a.dry_run else 'LIVE'}\n")

    if a.dry_run:
        print("dry run: plan only, nothing launched")
        return 0

    t0 = time.time()
    log = [{"kind": "meta", "config": cfg, "hardware": hw,
            "authorized_budget": authorized, "started": time.strftime("%Y-%m-%d %H:%M:%S")}]
    rng = random.Random(cfg["seeds"][0])
    run_verifier(cfg, hw, rng, log)
    print(f"verifier records: {sum(1 for r in log if r['kind']=='verifier')} "
          f"({time.time()-t0:.1f}s)")
    try:
        run_predictor(cfg, hw, log)
    except Exception as exc:
        log.append({"kind": "predictor_error", "error": f"{type(exc).__name__}: {exc}"})
        print(f"predictor arm failed: {type(exc).__name__}: {exc}")
    elapsed = time.time() - t0
    core_hours = elapsed / 3600.0
    log.append({"kind": "final", "elapsed_s": elapsed, "core_hours": core_hours,
                "within_budget": core_hours <= authorized["cpu_core_hours"],
                "wall_within_budget": elapsed / 3600.0 <= authorized["wall_clock_hours"]})

    out = pilot.get("out_dir", "runs/pilot")
    os.makedirs(out, exist_ok=True)
    path = os.path.join(out, f"pilot-{time.strftime('%Y%m%d-%H%M%S')}.json")
    with open(path, "w") as fh:
        json.dump(log, fh, indent=1, default=str)

    nv = sum(1 for r in log if r["kind"] == "verifier")
    npd = sum(1 for r in log if r["kind"] == "predictor")
    nvoid = sum(1 for r in log if r["kind"] == "void")
    nto = sum(1 for r in log if r.get("timeout"))
    nerr = sum(1 for r in log if r["kind"] in ("error", "predictor_error"))
    print(f"predictor records: {npd}")
    print(f"void (excluded)  : {nvoid}   timeouts: {nto}   errors: {nerr}")
    print(f"elapsed          : {elapsed:.1f}s ({core_hours:.4f} core-hours) "
          f"of {authorized['cpu_core_hours']} authorised")
    print(f"log              : {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
