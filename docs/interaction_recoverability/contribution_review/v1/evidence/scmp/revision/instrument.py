"""Component-level instrumentation for the verifier.

Protocol v1 recorded one terminal `wall_s` per instance, which is a STOPPING time
and for a timed-out record equals the budget. That cannot answer where the time
goes, and it cannot give a first-resolution time. This module wraps the seven
components the plan names and records, per instance:

  chart_compilation   oracle construction and any per-family precomputation
  bound_propagation   affine envelope propagation through the network
  support_query       oracle support/count/feasibility calls
  cache_access        lookups, split into hits and misses, with the work each did
  policy_inference    learned scorer forward passes
  incumbent_search    feasible-point evaluation
  branching           split construction and child bookkeeping

It also records a (t, L, U, status) TRACE at every bound update, which is what
protocol v2 needs to derive first-resolution time and an anytime gap AUC. The
wrappers count and time; they never change a returned value, and a test asserts
that instrumented and uninstrumented runs agree exactly.
"""
from __future__ import annotations

import time
from contextlib import contextmanager

COMPONENTS = ("chart_compilation", "bound_propagation", "support_query",
              "cache_access", "policy_inference", "incumbent_search", "branching")


class Meter:
    """Wall time, call count and 'actual work' per component, plus a bound trace."""

    def __init__(self, tolerance=1e-9):
        self.t0 = time.perf_counter()
        self.time = {c: 0.0 for c in COMPONENTS}
        self.calls = {c: 0 for c in COMPONENTS}
        self.work = {c: 0 for c in COMPONENTS}
        self.cache = {"hits": 0, "misses": 0, "hit_time": 0.0, "miss_time": 0.0,
                      "work_on_hit": 0, "work_on_miss": 0}
        self.trace = []                  # (t, lower, upper, status)
        self.tolerance = tolerance
        self.first_resolution_s = None
        self.first_incumbent_s = None

    # ------------------------------------------------------------------ timing
    @contextmanager
    def timed(self, component, work=0):
        assert component in COMPONENTS, component
        s = time.perf_counter()
        try:
            yield
        finally:
            d = time.perf_counter() - s
            self.time[component] += d
            self.calls[component] += 1
            self.work[component] += int(work)

    @contextmanager
    def cache_lookup(self):
        """Charge a cache access and record whether it did real work.

        A hit that avoids an oracle call is cheap; a miss that triggers one is
        not. Protocol v1 charged both alike, which made the budget counter read
        80 passes against a budget of 8.
        """
        s = time.perf_counter()
        box = {"hit": False, "work": 0}
        try:
            yield box
        finally:
            d = time.perf_counter() - s
            self.time["cache_access"] += d
            self.calls["cache_access"] += 1
            if box["hit"]:
                self.cache["hits"] += 1
                self.cache["hit_time"] += d
                self.cache["work_on_hit"] += int(box["work"])
            else:
                self.cache["misses"] += 1
                self.cache["miss_time"] += d
                self.cache["work_on_miss"] += int(box["work"])

    # -------------------------------------------------------------- bound trace
    def note_bound(self, lower, upper, status="open"):
        """Record a bound update. This is what first-resolution time is read from."""
        t = time.perf_counter() - self.t0
        self.trace.append({"t": t, "lower": float(lower), "upper": float(upper),
                           "status": status})
        if self.first_incumbent_s is None and lower > float("-inf"):
            self.first_incumbent_s = t
        if (self.first_resolution_s is None
                and upper - lower <= self.tolerance
                and lower > float("-inf")):
            self.first_resolution_s = t

    # ------------------------------------------------------------------ metrics
    def elapsed(self):
        return time.perf_counter() - self.t0

    def anytime_auc(self, horizon, scale):
        """(1/T) * integral_0^T min(1, gap(t)/scale) dt, from the step trace.

        `scale` is a common per-instance value shared by every arm, so the AUCs
        are comparable. Before the first valid bound the integrand is 1 -- the
        arm has produced nothing, and substituting zero there would reward
        silence. A negative gap is a defect: it is recorded and clamped at 0
        only for the integral, never silently repaired in the record.
        """
        if scale <= 0 or horizon <= 0:
            return float("nan"), 0
        pts = [p for p in self.trace if p["t"] <= horizon]
        area, prev_t, prev_val, negatives = 0.0, 0.0, 1.0, 0
        for p in pts:
            area += prev_val * (p["t"] - prev_t)
            gap = p["upper"] - p["lower"]
            if gap < 0:
                negatives += 1
            prev_val = min(1.0, max(0.0, gap) / scale) if gap == gap else 1.0
            prev_t = p["t"]
        area += prev_val * (horizon - prev_t)
        return area / horizon, negatives

    def summary(self):
        tot = sum(self.time.values())
        return {
            "component_time_s": {k: round(v, 6) for k, v in self.time.items()},
            "component_share": {k: (round(v / tot, 4) if tot > 0 else 0.0)
                                for k, v in self.time.items()},
            "component_calls": dict(self.calls),
            "component_work": dict(self.work),
            "cache": dict(self.cache),
            "cache_hit_rate": (self.cache["hits"]
                               / max(self.cache["hits"] + self.cache["misses"], 1)),
            "instrumented_total_s": round(tot, 6),
            "wall_s": round(self.elapsed(), 6),
            "first_incumbent_s": self.first_incumbent_s,
            "first_resolution_s": self.first_resolution_s,
            "trace_points": len(self.trace),
        }
