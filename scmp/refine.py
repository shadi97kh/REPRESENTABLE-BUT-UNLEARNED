"""Branch refinement over the feasible family, with a replayable certificate.

SPLITS COVER. Splitting node `S` on candidate edge `e` produces `S & {z_e=1}` and
`S & {z_e=0}`. Every member of `S` satisfies exactly one, so the children cover the
parent with no gap and no overlap. Nothing can be lost by splitting.

CHILD CAP = min(parent, new). A child's family is a subset of its parent's, so the
parent's bound is already valid for it. Taking the minimum means a child bound can only
improve on what was inherited, which is what makes the global bound monotone.

PRUNING. A node is pruned only when a VALID upper bound on it cannot beat the incumbent.
The incumbent is always a real model evaluation at a feasible point, so pruning discards
only subtrees that provably contain nothing better.

ELIMINATION. A child is removed only when the oracle PROVES its family empty. Ranking a
branch low never removes it: the policy chooses the order of work, never its extent.

FAIRNESS. Pure best-first can starve a node whose bound is mediocre but which is never
the best. Every `fairness_period` expansions the queue is served FIFO instead, so every
open non-singleton is expanded eventually.

NUMERICAL STATUS is inherited from the bounds it consumes. If those are
`float_unverified` the result is `float_unverified`: refinement cannot manufacture
rigour that the arithmetic beneath it does not have.
"""
from __future__ import annotations

import heapq
import itertools
import json
import time
from dataclasses import dataclass, field

import torch

from .bounds import certified_upper_bound, _model_hash
from .numerics import FLOAT_UNVERIFIED, combine_status

PROVED = "proved"
UNRESOLVED = "unresolved"
INFEASIBLE = "infeasible"


@dataclass
class Node:
    id: int
    parent: int | None
    forced: tuple
    forbidden: tuple
    cap: float
    singleton: bool = False


@dataclass
class RefineResult:
    status: str
    global_upper: float
    incumbent: float
    incumbent_witness: tuple
    gap: float
    numerical_status: str
    records: list = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    def certificate(self) -> dict:
        return {"schema": "scmp.refine/1", "status": self.status,
                "global_upper": self.global_upper, "incumbent": self.incumbent,
                "incumbent_witness": [list(e) for e in self.incumbent_witness],
                "gap": self.gap, "numerical_status": self.numerical_status,
                "records": self.records, "stats": self.stats}


def _undecided(cands, forced, forbidden):
    f, b = set(forced), set(forbidden)
    return [e for e in cands if e not in f and e not in b]


def refine(model, X, oracle, base=None, forced=(), forbidden=(), *,
           max_nodes: int = 200, time_limit: float = 30.0, tol: float = 1e-9,
           fairness_period: int = 4, reference_mode: bool = True,
           bound_fn=None, split_policy=None, seed: int = 0) -> RefineResult:
    from .model import backbone_adjacency
    if base is None:
        base = backbone_adjacency(model.n_nodes, dtype=model.dtype)
    bound_fn = bound_fn or (lambda f, b: certified_upper_bound(
        model, X, oracle, base, f, b, mode="combined"))

    t0 = time.time()
    records = []
    model_hash = _model_hash(model)
    root_feas = oracle.feasibility(forced, forbidden)
    records.append({"op": "root", "model_hash": model_hash,
                    "family_hash": root_feas.family_hash,
                    "forced": [list(e) for e in forced],
                    "forbidden": [list(e) for e in forbidden]})
    if not root_feas.feasible:
        return RefineResult(INFEASIBLE, float("-inf"), float("-inf"), (), 0.0,
                            FLOAT_UNVERIFIED, records,
                            {"reason": root_feas.reason})

    incumbent, incumbent_witness = float("-inf"), ()
    statuses = []

    def evaluate(witness):
        with torch.no_grad():
            return float(model(X, model.indicator(witness), base))

    def offer(witness, node_id, why):
        """Feasible incumbents only: a real model evaluation at a feasible point."""
        nonlocal incumbent, incumbent_witness
        val = evaluate(witness)
        if val > incumbent:
            incumbent, incumbent_witness = val, tuple(sorted(witness))
            records.append({"op": "incumbent", "node": node_id, "value": val,
                            "witness": [list(e) for e in incumbent_witness],
                            "why": why})
        return val

    counter = itertools.count()
    root_b = bound_fn(tuple(forced), tuple(forbidden))
    statuses.append(root_b.numerical_status)
    offer(root_b.incumbent_witness, 0, "root bound witness")
    root = Node(0, None, tuple(forced), tuple(forbidden), float(root_b.upper))
    records.append({"op": "bound", "node": 0, "parent": None,
                    "raw": float(root_b.upper), "parent_cap": None,
                    "cap": root.cap})

    heap = [(-root.cap, next(counter), root)]
    fifo = [root]
    open_nodes = {0: root}
    next_id = 1
    expansions = 0
    global_upper = root.cap
    history = [global_upper]
    stats = {"expanded": 0, "pruned": 0, "eliminated": 0, "singletons": 0,
             "bounds": 1, "fairness_pops": 0}

    while heap or fifo:
        if time.time() - t0 > time_limit or expansions >= max_nodes:
            break
        use_fifo = fairness_period and expansions and expansions % fairness_period == 0
        node = None
        if use_fifo:
            while fifo:
                cand = fifo.pop(0)
                if cand.id in open_nodes:
                    node = cand
                    stats["fairness_pops"] += 1
                    break
        if node is None:
            while heap:
                _, _, cand = heapq.heappop(heap)
                if cand.id in open_nodes:
                    node = cand
                    break
        if node is None:
            break
        open_nodes.pop(node.id, None)

        if node.cap <= incumbent + tol:
            records.append({"op": "prune", "node": node.id, "cap": node.cap,
                            "incumbent": incumbent})
            stats["pruned"] += 1
            global_upper = max(incumbent, max((n.cap for n in open_nodes.values()),
                                              default=float("-inf")))
            history.append(global_upper)
            continue

        undecided = _undecided(model.candidates, node.forced, node.forbidden)
        cnt = oracle.count(node.forced, node.forbidden).value
        if reference_mode and cnt == 1:
            w = oracle.support({e: 0.0 for e in model.candidates},
                               node.forced, node.forbidden).witness
            val = offer(w, node.id, "exact singleton evaluation")
            records.append({"op": "singleton", "node": node.id, "count": cnt,
                            "value": val, "witness": [list(e) for e in w]})
            stats["singletons"] += 1
            global_upper = max(incumbent, max((n.cap for n in open_nodes.values()),
                                              default=float("-inf")))
            history.append(global_upper)
            continue
        if not undecided:
            records.append({"op": "leaf", "node": node.id, "count": cnt})
            global_upper = max(incumbent, max((n.cap for n in open_nodes.values()),
                                              default=float("-inf")))
            history.append(global_upper)
            continue

        edge = (split_policy(undecided, node) if split_policy else undecided[0])
        expansions += 1
        stats["expanded"] += 1
        child_ids = []
        for kind in ("forced", "forbidden"):
            cf = tuple(sorted(set(node.forced) | {edge})) if kind == "forced" else node.forced
            cb = node.forbidden if kind == "forced" else tuple(sorted(set(node.forbidden) | {edge}))
            feas = oracle.feasibility(cf, cb)
            cid = next_id
            next_id += 1
            child_ids.append(cid)
            if not feas.feasible:
                records.append({"op": "eliminate", "node": cid, "parent": node.id,
                                "side": kind, "edge": list(edge),
                                "reason": feas.reason})
                stats["eliminated"] += 1
                continue
            rb = bound_fn(cf, cb)
            statuses.append(rb.numerical_status)
            stats["bounds"] += 1
            cap = min(node.cap, float(rb.upper))          # child cap = min(parent, new)
            records.append({"op": "bound", "node": cid, "parent": node.id,
                            "side": kind, "edge": list(edge),
                            "forced": [list(e) for e in cf],
                            "forbidden": [list(e) for e in cb],
                            "raw": float(rb.upper), "parent_cap": node.cap,
                            "cap": cap})
            offer(rb.incumbent_witness, cid, "child bound witness")
            if cap <= incumbent + tol:
                records.append({"op": "prune", "node": cid, "cap": cap,
                                "incumbent": incumbent})
                stats["pruned"] += 1
                continue
            child = Node(cid, node.id, cf, cb, cap)
            open_nodes[cid] = child
            heapq.heappush(heap, (-cap, next(counter), child))
            fifo.append(child)
        records.append({"op": "split", "node": node.id, "edge": list(edge),
                        "children": child_ids})

        global_upper = max(incumbent, max((n.cap for n in open_nodes.values()),
                                          default=float("-inf")))
        history.append(global_upper)

    resolved = not open_nodes and not (time.time() - t0 > time_limit
                                       or expansions >= max_nodes)
    status = PROVED if resolved else UNRESOLVED
    if resolved:
        global_upper = incumbent
        history.append(global_upper)
    stats.update({"nodes": next_id, "open": len(open_nodes),
                  "elapsed_s": time.time() - t0, "history": history})
    records.append({"op": "final", "status": status, "global_upper": global_upper,
                    "incumbent": incumbent})
    return RefineResult(status, float(global_upper), float(incumbent),
                        incumbent_witness, float(global_upper - incumbent),
                        combine_status(*statuses) if statuses else FLOAT_UNVERIFIED,
                        records, stats)
