"""A development-only difficulty ladder with a FROZEN generation rule.

The v1 pilot families were closed to gap zero by the plain box baseline on 99.1%
of instances, so G1 and G2 could not discriminate the arms at all. This module
builds a graded set instead. Two rules govern it:

  * the rule is declared and hashed BEFORE any evaluation case is drawn, and the
    hash is written into every record. Cases are never selected because an arm
    wins on them, and never dropped for being easy.
  * the easy v1 cases are PRESERVED as rung 0 and always appear in aggregate
    reporting. Reporting only the hard stratum would be selection by another name.

Difficulty is graded along the five axes the plan names -- candidate-edge count,
family compatibility (how much structure the constraints force), supported
depth/width, unstable ReLU count and decision margin -- so a rung is a
declared point in that space, not a post-hoc label.
"""
from __future__ import annotations

import hashlib
import json

# The frozen rule. Editing this dict changes RULE_HASH, which invalidates every
# record made under the old rule rather than silently mixing generations.
LADDER = {
    "version": "v2.0",
    "rungs": [
        {"rung": 0, "name": "v1_easy_preserved",
         "families": ["path_matching:6", "layered_dag:3x3", "rna:GGAAACC"],
         "d_hid": 4, "perturb": 1.0, "gamma": 1.0,
         "note": "the protocol-v1 instances, kept verbatim so aggregate "
                 "reporting still contains what v1 measured"},
        {"rung": 1, "name": "wider_candidate_set",
         "families": ["path_matching:8", "layered_dag:3x4", "rna:GCGCAAAAGCGC"],
         "d_hid": 6, "perturb": 1.5, "gamma": 1.0,
         "axis": "candidate-edge count"},
        {"rung": 2, "name": "looser_compatibility",
         "families": ["path_matching:10", "layered_dag:4x4", "rna:GGGGAAAACCCC"],
         "d_hid": 8, "perturb": 2.0, "gamma": 1.0,
         "axis": "family compatibility: more members per candidate edge"},
        {"rung": 3, "name": "wider_network_more_unstable_relus",
         "families": ["path_matching:10", "layered_dag:4x5", "rna:GGCGCAAAAGCGCC"],
         "d_hid": 12, "perturb": 3.0, "gamma": 1.0,
         "axis": "supported width and unstable-ReLU count"},
    ],
    "seeds_per_rung": 3,
    "decision_margin_grid": [0.02, 0.10, 0.50],
    "margin_definition":
        "tau is placed at incumbent + margin * |incumbent|, so a decision is "
        "'is max f over the family below tau'. Small margins are near-boundary "
        "and hard; the grid is declared here, not chosen after seeing results.",
    "selection_policy":
        "every generated case is evaluated and reported. No case is added, "
        "dropped or reweighted after any arm has been run on it.",
}


def rule_hash() -> str:
    return hashlib.sha256(
        json.dumps(LADDER, sort_keys=True).encode()).hexdigest()[:16]


RULE_HASH = rule_hash()


def _oracle_for(tag):
    from ..oracles.layered_dag import LayeredDAGPathOracle
    from ..oracles.path_matching import PathMatchingOracle
    from ..oracles.rna_noncrossing import RNANonCrossingOracle
    kind, arg = tag.split(":", 1)
    if kind == "path_matching":
        return PathMatchingOracle(int(arg)), int(arg)
    if kind == "layered_dag":
        # The oracle's ground set is (layer, u, v) triples, but CertifiedModel
        # indexes node PAIRS. EdgeRelabeledOracle is the existing bridge and is
        # reused here exactly as run_pilot.build_family uses it, so the two
        # stages present the same family under the same names.
        from ..oracles.adapters import EdgeRelabeledOracle
        L, W = (int(x) for x in arg.split("x"))
        sizes = tuple([W] * L)
        edges = [(l, u, v) for l in range(L - 1)
                 for u in range(W) for v in range(W)]
        base = LayeredDAGPathOracle(sizes, edges)
        flat, i = {}, 0
        for li in range(L):
            for u in range(W):
                flat[(li, u)] = i
                i += 1
        mapping = {e: tuple(sorted((flat[(e[0], e[1])], flat[(e[0] + 1, e[2])])))
                   for e in base.ground_set()}
        return EdgeRelabeledOracle(base, mapping), i
    if kind == "rna":
        return RNANonCrossingOracle(arg, 3, True), len(arg)
    raise ValueError(tag)


def generate(rungs=None, seeds_per_rung=None, margins=None):
    """Yield every case the frozen rule specifies. Deterministic and total."""
    out = []
    for spec in LADDER["rungs"]:
        if rungs is not None and spec["rung"] not in rungs:
            continue
        for tag in spec["families"]:
            oracle, n_nodes = _oracle_for(tag)
            for seed in range(seeds_per_rung or LADDER["seeds_per_rung"]):
                for margin in (margins or LADDER["decision_margin_grid"]):
                    out.append({
                        "rung": spec["rung"], "rung_name": spec["name"],
                        "family_tag": tag, "n_nodes": n_nodes,
                        "n_candidates": len(oracle.ground_set()),
                        "d_hid": spec["d_hid"], "perturb": spec["perturb"],
                        "gamma": spec["gamma"], "seed": seed, "margin": margin,
                        "rule_hash": RULE_HASH})
    return out


def build_case(case, d_in=4):
    """Instantiate the model, features and oracle for one declared case."""
    import torch

    from ..model import CertifiedModel, backbone_adjacency
    oracle, n_nodes = _oracle_for(case["family_tag"])
    edges = oracle.ground_set()
    m = CertifiedModel(d_in, edges, n_nodes, d_hid=case["d_hid"],
                       gamma=case["gamma"], seed=case["seed"])
    g = torch.Generator().manual_seed(1000 + case["seed"])
    with torch.no_grad():
        for p in m.parameters():
            p.add_(torch.randn(p.shape, generator=g, dtype=p.dtype)
                   * case["perturb"])
    X = torch.rand(n_nodes, d_in, dtype=torch.float64,
                   generator=torch.Generator().manual_seed(case["seed"]))
    return m, X, oracle, backbone_adjacency(n_nodes)
