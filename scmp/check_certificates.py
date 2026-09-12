"""Independent checker for refinement certificates.

This module does NOT import `scmp.refine`. It re-derives every claim from the record
stream and the model, so a bug in the search cannot excuse itself. Where the family is
small enough it also enumerates it, which turns "the recorded bound is consistent" into
"the recorded bound is correct".

Run:  python -m scmp.check_certificates --run runs/tiny_complete
"""
from __future__ import annotations

import argparse
import glob
import json
import os
import sys

import torch

from .bounds import _model_hash
from .model import CertifiedModel, backbone_adjacency
from .oracles import RNANonCrossingOracle
from .oracles.bruteforce import enumerate_rna

TOL = 1e-9


def _t(pairs):
    return tuple(sorted(tuple(e) for e in pairs))


class Checker:
    def __init__(self, cert, model, X, oracle, base, members=None, seq=None):
        self.cert, self.model, self.X = cert, model, X
        self.oracle, self.base = oracle, base
        self.members = members
        self.seq = seq
        self.failures = []
        self.checks = 0

    def fail(self, kind, detail):
        self.failures.append({"kind": kind, **detail})

    def ok(self):
        self.checks += 1

    # -- individual obligations -------------------------------------------
    def check_identity(self):
        root = next(r for r in self.cert["records"] if r["op"] == "root")
        self.ok()
        if root["model_hash"] != _model_hash(self.model):
            self.fail("stale_model_hash", {"cert": root["model_hash"],
                                           "actual": _model_hash(self.model)})
        feas = self.oracle.feasibility(_t(root["forced"]), _t(root["forbidden"]))
        self.ok()
        if root["family_hash"] != feas.family_hash:
            self.fail("stale_family_hash", {"cert": root["family_hash"],
                                            "actual": feas.family_hash})

    def check_coverage(self):
        """Children of a split must be parent+{e=1} and parent+{e=0}: exact cover."""
        masks = {}
        for r in self.cert["records"]:
            if r["op"] == "root":
                masks[0] = (_t(r["forced"]), _t(r["forbidden"]))
            elif r["op"] == "bound" and r.get("parent") is not None:
                masks[r["node"]] = (_t(r["forced"]), _t(r["forbidden"]))
        by_parent = {}
        for r in self.cert["records"]:
            if r["op"] == "bound" and r.get("parent") is not None:
                by_parent.setdefault(r["parent"], []).append(r)
            if r["op"] == "eliminate":
                by_parent.setdefault(r["parent"], []).append(r)
        for split in [r for r in self.cert["records"] if r["op"] == "split"]:
            self.ok()
            pid, edge = split["node"], tuple(split["edge"])
            if pid not in masks:
                self.fail("split_of_unknown_node", {"node": pid})
                continue
            pf, pb = masks[pid]
            want = {("forced", _t(set(pf) | {edge}), pb),
                    ("forbidden", pf, _t(set(pb) | {edge}))}
            got = set()
            for ch in by_parent.get(pid, []):
                if ch["op"] == "bound":
                    got.add((ch["side"], _t(ch["forced"]), _t(ch["forbidden"])))
                else:
                    side = ch["side"]
                    cf = _t(set(pf) | {edge}) if side == "forced" else pf
                    cb = pb if side == "forced" else _t(set(pb) | {edge})
                    got.add((side, cf, cb))
            if got != want:
                self.fail("split_does_not_cover",
                          {"node": pid, "edge": list(edge),
                           "expected": sorted(map(str, want)),
                           "found": sorted(map(str, got))})

    def check_caps(self):
        caps = {}
        for r in self.cert["records"]:
            if r["op"] != "bound":
                continue
            self.ok()
            if r.get("parent") is None:
                caps[r["node"]] = r["cap"]
                if r["cap"] != r["raw"]:
                    self.fail("root_cap_mismatch", {"node": r["node"]})
                continue
            want = min(r["parent_cap"], r["raw"])
            if abs(r["cap"] - want) > 0:
                self.fail("cap_not_min_of_parent_and_new",
                          {"node": r["node"], "cap": r["cap"], "expected": want})
            if r["parent"] in caps and r["cap"] > caps[r["parent"]] + TOL:
                self.fail("cap_exceeds_parent", {"node": r["node"]})
            caps[r["node"]] = r["cap"]

    def check_eliminations(self):
        """Re-prove every elimination against the oracle."""
        masks = {}
        for r in self.cert["records"]:
            if r["op"] == "root":
                masks[0] = (_t(r["forced"]), _t(r["forbidden"]))
            elif r["op"] == "bound" and r.get("parent") is not None:
                masks[r["node"]] = (_t(r["forced"]), _t(r["forbidden"]))
        for r in [x for x in self.cert["records"] if x["op"] == "eliminate"]:
            self.ok()
            pf, pb = masks.get(r["parent"], ((), ()))
            edge = tuple(r["edge"])
            cf = _t(set(pf) | {edge}) if r["side"] == "forced" else pf
            cb = pb if r["side"] == "forced" else _t(set(pb) | {edge})
            if self.oracle.feasibility(cf, cb).feasible:
                self.fail("eliminated_a_feasible_branch",
                          {"node": r["node"], "forced": list(cf),
                           "forbidden": list(cb)})

    def check_incumbents(self):
        for r in [x for x in self.cert["records"] if x["op"] in ("incumbent",
                                                                 "singleton")]:
            self.ok()
            w = _t(r["witness"])
            if self.members is not None and w not in {_t(s) for s in self.members}:
                self.fail("incumbent_infeasible", {"witness": list(w)})
            with torch.no_grad():
                direct = float(self.model(self.X, self.model.indicator(w), self.base))
            if abs(direct - r["value"]) > 0:
                self.fail("incumbent_not_a_model_evaluation",
                          {"witness": list(w), "recorded": r["value"],
                           "recomputed": direct})

    def check_prunes(self):
        for r in [x for x in self.cert["records"] if x["op"] == "prune"]:
            self.ok()
            if r["cap"] > r["incumbent"] + TOL:
                self.fail("pruned_a_node_that_could_beat_the_incumbent",
                          {"node": r["node"], "cap": r["cap"],
                           "incumbent": r["incumbent"]})

    def check_monotone_global(self):
        hist = self.cert.get("stats", {}).get("history", [])
        for a, b in zip(hist, hist[1:]):
            self.ok()
            if b > a + TOL:
                self.fail("global_upper_increased", {"from": a, "to": b})

    def check_against_enumeration(self):
        """Only possible on tiny families, and the strongest check available."""
        if self.members is None:
            return
        with torch.no_grad():
            true_max = max(float(self.model(self.X, self.model.indicator(s), self.base))
                           for s in self.members)
        self.ok()
        if self.cert["global_upper"] < true_max - TOL:
            self.fail("global_upper_below_true_max",
                      {"global_upper": self.cert["global_upper"],
                       "true_max": true_max})
        self.ok()
        if self.cert["incumbent"] > true_max + TOL:
            self.fail("incumbent_above_true_max",
                      {"incumbent": self.cert["incumbent"], "true_max": true_max})
        if self.cert["status"] == "proved":
            self.ok()
            if abs(self.cert["global_upper"] - true_max) > 1e-9:
                self.fail("claimed_proved_but_not_exact",
                          {"global_upper": self.cert["global_upper"],
                           "true_max": true_max})

    def check_status_honesty(self):
        self.ok()
        if self.cert["numerical_status"] not in ("float_unverified", "float_enclosed",
                                                 "exact_rational"):
            self.fail("unknown_numerical_status",
                      {"status": self.cert["numerical_status"]})
        if self.cert["status"] == "unresolved":
            self.ok()
            if self.cert["gap"] < -TOL:
                self.fail("unresolved_with_negative_gap", {"gap": self.cert["gap"]})

    def run(self):
        self.check_identity()
        self.check_coverage()
        self.check_caps()
        self.check_eliminations()
        self.check_incumbents()
        self.check_prunes()
        self.check_monotone_global()
        self.check_against_enumeration()
        self.check_status_honesty()
        return {"checks": self.checks, "failures": self.failures,
                "accepted": not self.failures}


def rebuild(meta):
    """Reconstruct the exact model the certificate was produced for."""
    seq = meta["seq"]
    o = RNANonCrossingOracle(seq, meta["min_loop"], meta["canonical_only"])
    m = CertifiedModel(4, o.ground_set(), len(seq), d_hid=meta["d_hid"],
                       gamma=meta["gamma"], seed=meta["init_seed"])
    if meta["perturb_scale"]:
        g = torch.Generator().manual_seed(meta["perturb_seed"])
        with torch.no_grad():
            for p in m.parameters():
                p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype)
                       * meta["perturb_scale"])
    X = torch.zeros(len(seq), 4, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
    return o, m, X, backbone_adjacency(len(seq))


def check_file(path):
    payload = json.load(open(path))
    meta, cert = payload["meta"], payload["certificate"]
    o, m, X, B = rebuild(meta)
    members = enumerate_rna(meta["seq"], meta["min_loop"], meta["canonical_only"])
    return Checker(cert, m, X, o, B, members, meta["seq"]).run()


def main(argv=None):
    ap = argparse.ArgumentParser(prog="python -m scmp.check_certificates")
    ap.add_argument("--run", required=True)
    a = ap.parse_args(argv)
    files = sorted(glob.glob(os.path.join(a.run, "*.json")))
    if not files:
        print(f"no certificates under {a.run}")
        return 1
    print(f"independent certificate check: {len(files)} file(s) under {a.run}")
    print("(this module never imports scmp.refine)\n")
    total_checks = 0
    bad = []
    for f in files:
        rep = check_file(f)
        total_checks += rep["checks"]
        mark = "accept" if rep["accepted"] else "REJECT"
        print(f"  {mark}  {os.path.basename(f):40s} {rep['checks']:5d} obligations")
        if not rep["accepted"]:
            bad.append((f, rep))
            for fl in rep["failures"][:5]:
                print(f"          {fl}")
    print(f"\nobligations discharged : {total_checks}")
    print(f"certificates accepted  : {len(files) - len(bad)}/{len(files)}")
    if bad:
        print("\nREJECTED certificates present")
        return 2
    print("all certificates accepted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
