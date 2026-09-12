"""Matchings on a path graph: the second non-RNA family.

FAMILY. Vertices 0..n-1, candidate edges (i, i+1). A member is a matching: no two edges
share a vertex, which for a path means no two consecutive edge indices.

WHY THIS FAMILY. It is the one C5 used to show that relaxation slack is unbounded --
matchings are downward closed but not union closed, so the smallest containing lattice
has an infeasible top element. It is edge-indexed, so it plugs into `CertifiedModel`
unchanged, and it shares no structure with non-crossing pairing beyond being a family of
edge subsets. If an assumption about RNA has leaked into the bound machinery, running the
same pipeline here is what exposes it.

RECURRENCE, leftmost-anchored and unambiguous over edge index i:

    W[i] = W[i+1]                  ... edge i unused
         + w_i * W[i+2]            ... edge i used, so edge i+1 cannot be
    W[i] = 1 for i >= n-1

Each matching has exactly one derivation, keyed on the smallest used index. Counts are
Fibonacci numbers, which `tests/oracles` can check against a closed form.
"""
from __future__ import annotations

from .base import (COUNTING, EXACT_INTEGER, FLOAT_DIAGNOSTIC, FLOAT_UNVERIFIED,
                   INFEASIBLE, LOG_SUM_EXP, MAX_PLUS, NEG_INF, PROVED, Feasibility,
                   OracleResult, Production, SupportOracle, family_hash, normalise_mask)


class PathMatchingOracle(SupportOracle):
    family_id = "path_matching"

    def __init__(self, n_nodes: int):
        self.n = int(n_nodes)
        self._ground = tuple((i, i + 1) for i in range(max(self.n - 1, 0)))
        self._ground_set = frozenset(self._ground)
        self._last_stats: dict = {}

    def ground_set(self):
        return self._ground

    def describe(self):
        return {"family": self.family_id, "n_nodes": self.n,
                "edges": [list(e) for e in self._ground]}

    def assumptions(self):
        return ("path graph on n nodes", "matching: no two edges share a vertex",
                "candidate edges are exactly the path edges")

    def _hash(self, forced, forbidden):
        d = self.describe()
        d["forced"] = sorted(forced)
        d["forbidden"] = sorted(forbidden)
        return family_hash(d)

    def feasibility(self, forced=(), forbidden=()) -> Feasibility:
        forced, forbidden = normalise_mask(forced, forbidden)
        h = self._hash(forced, forbidden)
        A = self.assumptions()

        def no(reason):
            return Feasibility(False, reason, INFEASIBLE, h, A)

        both = forced & forbidden
        if both:
            return no(f"edge(s) both forced and forbidden: {sorted(both)}")
        for e in sorted(forced):
            if e not in self._ground_set:
                return no(f"forced edge {e} is not a path edge")
        fl = sorted(forced)
        for a in range(len(fl)):
            for b in range(a + 1, len(fl)):
                if set(fl[a]) & set(fl[b]):
                    return no(f"forced edges {fl[a]} and {fl[b]} share a vertex")
        c = self._run(COUNTING, {}, forced, forbidden)
        if c["value"] == 0:
            return no("no matching satisfies the mask (recurrence returned count 0)")
        return Feasibility(True, f"{c['value']} matching(s) satisfy the mask", PROVED,
                           h, A)

    def _run(self, semiring, coefficients, forced, forbidden, want_backtrace=False):
        m = len(self._ground)
        zero, one = semiring.zero, semiring.one
        add, mul = semiring.add, semiring.mul
        forced_idx = {self._ground.index(e) for e in forced if e in self._ground_set}
        forbid_idx = {self._ground.index(e) for e in forbidden if e in self._ground_set}

        W = [zero] * (m + 2)
        choice = [None] * (m + 2) if want_backtrace else None
        W[m] = one
        W[m + 1] = one
        chart_values = m + 2
        productions = 0

        for i in range(m - 1, -1, -1):
            val = zero
            best = None
            if i not in forced_idx:                      # production: edge i unused
                productions += 1
                val = W[i + 1]
                if want_backtrace and val != zero:
                    best = Production("skip", (i, i), None)
            if i not in forbid_idx:                      # production: edge i used
                if (i + 1) not in forced_idx:            # neighbour cannot also be forced
                    productions += 1
                    inner = W[i + 2]
                    if inner != zero:
                        w = semiring.weight(coefficients.get(self._ground[i], 0))
                        cand = mul(w, inner)
                        if want_backtrace and (val == zero or cand > val):
                            best = Production("use", (i, i), self._ground[i])
                        val = add(val, cand)
            W[i] = val
            if want_backtrace:
                choice[i] = best

        self._last_stats = {"chart_values": chart_values,
                            "materialized_productions": productions,
                            "n": self.n, "ground_set_size": m}
        return {"value": W[0] if m else one, "W": W, "choice": choice,
                "stats": dict(self._last_stats)}

    def _backtrace(self, choice):
        out, i, m = [], 0, len(self._ground)
        while i < m:
            c = choice[i]
            if c is None:
                break
            out.append(c)
            i += 1 if c.kind == "skip" else 2
        return out

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
        prods = self._backtrace(r["choice"]) if self._ground else []
        pairs = tuple(sorted(p.edge for p in prods if p.edge is not None))
        st = dict(r["stats"]); st["productions"] = tuple(prods)
        return OracleResult(r["value"], pairs, status, PROVED, h, "exact",
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

    def stats(self):
        return dict(self._last_stats)
