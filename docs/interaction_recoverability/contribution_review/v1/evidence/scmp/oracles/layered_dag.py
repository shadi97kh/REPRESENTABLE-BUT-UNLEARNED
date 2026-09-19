"""Source-to-sink paths in a layered DAG, through the same oracle interface.

A second, structurally unrelated family exists so that the interface in `base.py` is
tested against something other than the one family it was written for. If an assumption
about RNA leaked into the abstraction, this is where it shows up.

FAMILY. Layers of sizes [n_0, ..., n_T]. An edge is (l, u, v): node u in layer l to node
v in layer l+1. A member is the edge set of one path from some layer-0 node to some
layer-T node, so every member has exactly T edges, one per layer transition.

RECURRENCE. V[0][v] = 1 for every layer-0 node; V[l+1][w] = sum over edges (l, v, w) of
V[l][v] * w(l,v,w). Unambiguous: a path is determined by its edge sequence and each edge
is visited once, in layer order.

P maps each production (l, v, w) to the edge (l, v, w) itself — the identity, which makes
this family a useful control for the production-to-edge tests.

COMPLEXITY. Chart values sum(n_l); materialised productions |E|. Both counted at run time.
"""
from __future__ import annotations

from .base import (COUNTING, EXACT_INTEGER, FLOAT_DIAGNOSTIC, FLOAT_UNVERIFIED,
                   INFEASIBLE, LOG_SUM_EXP, MAX_PLUS, NEG_INF, PROVED, Feasibility,
                   OracleResult, Production, SupportOracle, family_hash, normalise_mask)


class LayeredDAGPathOracle(SupportOracle):
    family_id = "layered_dag_path"

    def __init__(self, layer_sizes, edges):
        self.layer_sizes = tuple(int(x) for x in layer_sizes)
        self.n_layers = len(self.layer_sizes)
        if self.n_layers < 2:
            raise ValueError("need at least two layers")
        seen = set()
        ground = []
        for e in edges:
            l, u, v = (int(x) for x in e)
            if not (0 <= l < self.n_layers - 1):
                raise ValueError(f"edge {e}: layer out of range")
            if not (0 <= u < self.layer_sizes[l] and 0 <= v < self.layer_sizes[l + 1]):
                raise ValueError(f"edge {e}: node out of range")
            t = (l, u, v)
            if t not in seen:
                seen.add(t); ground.append(t)
        self._ground = tuple(sorted(ground))
        self._ground_set = frozenset(self._ground)
        self._by_layer = {l: [] for l in range(self.n_layers - 1)}
        for (l, u, v) in self._ground:
            self._by_layer[l].append((l, u, v))
        self._last_stats: dict = {}

    def ground_set(self): return self._ground

    def describe(self):
        return {"family": self.family_id, "layer_sizes": list(self.layer_sizes),
                "edges": [list(e) for e in self._ground]}

    def assumptions(self):
        return ("paths run from a layer-0 node to a final-layer node",
                "exactly one edge per layer transition",
                "edges only between consecutive layers")

    def _hash(self, forced, forbidden):
        d = self.describe()
        d["forced"] = sorted(forced); d["forbidden"] = sorted(forbidden)
        return family_hash(d)

    def feasibility(self, forced=(), forbidden=()) -> Feasibility:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        A = self.assumptions()

        def no(reason): return Feasibility(False, reason, INFEASIBLE, h, A)

        both = forced & forbidden
        if both:
            return no(f"edge(s) both forced and forbidden: {sorted(both)}")
        for e in sorted(forced):
            if e not in self._ground_set:
                return no(f"forced edge {e} is not in the ground set")
        per_layer = {}
        for e in sorted(forced):
            l = e[0]
            if l in per_layer:
                return no(f"forced edges {per_layer[l]} and {e} share layer {l}")
            per_layer[l] = e
        # chaining is a genuine constraint: consecutive forced edges must meet
        for l in sorted(per_layer):
            if l + 1 in per_layer and per_layer[l][2] != per_layer[l + 1][1]:
                return no(f"forced edges {per_layer[l]} and {per_layer[l+1]} do not meet")
        c = self._run(COUNTING, {}, forced, forbidden)
        if c["value"] == 0:
            return no("no path satisfies the mask (recurrence returned count 0)")
        return Feasibility(True, f"{c['value']} path(s) satisfy the mask", PROVED, h, A)

    def _run(self, semiring, coefficients, forced, forbidden, want_backtrace=False):
        zero, one = semiring.zero, semiring.one
        add, mul = semiring.add, semiring.mul
        # Every forced edge must appear in the member. A path uses exactly one edge per
        # layer transition, so two forced edges sharing a layer make the family empty --
        # a conjunction, not a disjunction. Getting this wrong silently returns the
        # count of paths using *either* one.
        forced_at = {}
        for e in forced:
            if forced_at.setdefault(e[0], e) != e:
                self._last_stats = {"chart_values": 0, "materialized_productions": 0,
                                    "n_layers": self.n_layers,
                                    "ground_set_size": len(self._ground)}
                return {"value": zero, "V": None, "back": None, "sink": None,
                        "stats": dict(self._last_stats)}

        V = [[zero] * s for s in self.layer_sizes]
        back = [[None] * s for s in self.layer_sizes] if want_backtrace else None
        for v in range(self.layer_sizes[0]):
            V[0][v] = one
        chart_values = sum(self.layer_sizes)
        productions = 0

        for l in range(self.n_layers - 1):
            req = forced_at.get(l)
            for (ll, u, v) in self._by_layer[l]:
                if req is not None and (ll, u, v) != req:
                    continue
                if (ll, u, v) in forbidden:
                    continue
                productions += 1
                if V[l][u] == zero:
                    continue
                w = semiring.weight(coefficients.get((ll, u, v), 0))
                cand = mul(V[l][u], w)
                if want_backtrace and (V[l + 1][v] == zero or cand > V[l + 1][v]):
                    back[l + 1][v] = (ll, u, v)
                V[l + 1][v] = add(V[l + 1][v], cand)

        total = zero
        best_sink = None
        for v in range(self.layer_sizes[-1]):
            if V[-1][v] == zero:
                continue
            if want_backtrace and (total == zero or V[-1][v] > total):
                best_sink = v
            total = add(total, V[-1][v])

        self._last_stats = {"chart_values": chart_values,
                            "materialized_productions": productions,
                            "n_layers": self.n_layers,
                            "ground_set_size": len(self._ground)}
        return {"value": total, "V": V, "back": back, "sink": best_sink,
                "stats": dict(self._last_stats)}

    def support(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        coefficients = dict(coefficients or {})
        exact = all(isinstance(v, int) for v in coefficients.values())
        status = EXACT_INTEGER if exact else FLOAT_UNVERIFIED
        r = self._run(MAX_PLUS, coefficients, forced, forbidden, want_backtrace=True)
        if r["value"] == NEG_INF:
            return OracleResult(NEG_INF, None, status, INFEASIBLE, h, "none",
                                self.assumptions(), r["stats"])
        prods, v, l = [], r["sink"], self.n_layers - 1
        while l > 0:
            e = r["back"][l][v]
            prods.append(Production("step", (l - 1, l), e))
            v = e[1]
            l -= 1
        prods.reverse()
        edges = tuple(p.edge for p in prods)
        st = dict(r["stats"]); st["productions"] = tuple(prods)
        return OracleResult(r["value"], edges, status, PROVED, h, "exact",
                            self.assumptions(), st)

    def count(self, forced=(), forbidden=()) -> OracleResult:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        r = self._run(COUNTING, {}, forced, forbidden)
        term = PROVED if r["value"] > 0 else INFEASIBLE
        return OracleResult(r["value"], None, EXACT_INTEGER, term, h, "exact",
                            self.assumptions(), r["stats"])

    def log_partition(self, coefficients, forced=(), forbidden=()) -> OracleResult:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        r = self._run(LOG_SUM_EXP, dict(coefficients or {}), forced, forbidden)
        term = PROVED if r["value"] != NEG_INF else INFEASIBLE
        return OracleResult(r["value"], None, FLOAT_DIAGNOSTIC, term, h, "exact",
                            self.assumptions(), r["stats"])

    def stats(self): return dict(self._last_stats)
