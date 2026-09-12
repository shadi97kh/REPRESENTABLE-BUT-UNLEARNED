"""Conditional message support.

THE IDEA. In the layer-2 message the term `z_e * h_u` contributes only when `z_e = 1`.
So the relevant bound on `h_u` is not its unconditional maximum but its maximum over
graphs that actually contain `e`:

    q_e  :=  max { d^T W_e h_u(z)  :  z in F, z_e = 1 }

and then `z_e * h_u(z) <= z_e * q_e` for every `z in F`: at `z_e = 1` the point lies in
the conditioned family so `h_u <= q_e`, and at `z_e = 0` both sides are zero. Sound, and
never worse than the unconditional envelope because we retain

    q_e_used = min(unconditional bound, conditional bound)

so a conditional pass that fails to tighten cannot loosen anything.

WHY IT BITES ON THIS FAMILY. Conditioning on `e = 1` forbids every edge that crosses `e`
or shares an endpoint with it. Those edges then cannot contribute to `h_u`, which is the
carrier-gating separation. Independently, the tighter interval can flip an unstable ReLU
to stable, which is the ReLU separation. Both are reproduced in
tests/conditional_support.

IGNORING AN EDGE. An edge is dropped from the message only after the oracle PROVES
`F and {z_e = 1}` empty. Absence of evidence is never treated as proof: an edge we simply
did not spend budget on keeps its unconditional bound.

RECURSION. Computing `q_e` needs `h_u` bounded under the condition, which means
re-running the earlier layers with the mask `forced + {e}` -- the conditional recursion
moves to earlier layers. Layer 1 is exactly affine, so conditioning there is exact.
All recursive work is bounded by an explicit query budget and every pass is cached.

NUMERICAL STATUS. float64 throughout, no directed rounding: `float_unverified`,
diagnostic. See docs/bounds_proof.md.
"""
from __future__ import annotations

import random
from dataclasses import dataclass, field

import numpy as np
import torch

from .bounds import (BOX, FLOAT_UNVERIFIED, INFEASIBLE, PROVED, Aff, ResidualBounder,
                     Trace, _split, _tightest_lower, _tightest_upper, _model_hash,
                     mccormick, relu_bounds)


# --------------------------------------------------------------------- concretisation
class FamilyConc:
    """Concretise an affine form over the FEASIBLE family, not the box.

    max_{z in F} (c^T z + k) = k + support_F(c). Strictly tighter than the box whenever
    the family excludes the box maximiser, which for a non-crossing family is most of
    the time.
    """

    def __init__(self, oracle, cands, forced=(), forbidden=(), counter=None):
        self.oracle, self.cands = oracle, tuple(cands)
        self.forced, self.forbidden = tuple(sorted(forced)), tuple(sorted(forbidden))
        self.counter = counter if counter is not None else {"oracle_calls": 0}
        self._memo = {}

    def _support(self, coef):
        key = coef.tobytes()
        hit = self._memo.get(key)
        if hit is not None:
            return hit
        d = {e: float(coef[k]) for k, e in enumerate(self.cands)}
        r = self.oracle.support(d, self.forced, self.forbidden)
        self.counter["oracle_calls"] += 1
        val = r.value
        self._memo[key] = val
        return val

    def hi(self, a: Aff) -> float:
        if not np.any(a.coef):
            return a.const
        return a.const + self._support(a.coef)

    def lo(self, a: Aff) -> float:
        if not np.any(a.coef):
            return a.const
        return a.const - self._support(-a.coef)


# ------------------------------------------------------------------------ cache keying
@dataclass(frozen=True)
class CacheKey:
    model_hash: str
    family_hash: str
    mask: tuple
    layer: str
    node: int
    direction: int
    condition: tuple

    @staticmethod
    def make(model_hash, family_hash, forced, forbidden, layer, node, direction, cond):
        return CacheKey(model_hash, family_hash,
                        (tuple(sorted(forced)), tuple(sorted(forbidden))),
                        layer, int(node), int(direction), tuple(sorted(cond)))


class BoundCache:
    """Keyed by model, family, mask, layer, node, direction and condition.

    Invalidation is by construction: the model hash covers every weight and gamma, and
    the family hash covers the sequence, the min-loop and canonicity settings and the
    mask. Any change to any of those yields a different key, so a stale entry cannot be
    read back. There is no time-based or manual invalidation to get wrong.
    """

    def __init__(self):
        self.store = {}
        self.hits = 0
        self.misses = 0

    def get(self, key):
        if key in self.store:
            self.hits += 1
            return self.store[key]
        self.misses += 1
        return None

    def put(self, key, value):
        self.store[key] = value
        return value

    def stats(self):
        return {"entries": len(self.store), "hits": self.hits, "misses": self.misses}


# ---------------------------------------------------------------------------- q bound
@dataclass
class QBound:
    """A certified bound on d^T W_e h_u over valid graphs with e = 1."""
    edge: tuple
    node: int
    direction: int
    value: float
    source: str                    # conditional | unconditional_fallback | impossible
    unconditional: float
    conditional: float | None
    proved_impossible: bool
    assumptions: tuple = ()
    numerical_status: str = FLOAT_UNVERIFIED


# ----------------------------------------------------------------------- edge scorers
def score_random(edges, ctx, rng):
    order = list(edges)
    rng.shuffle(order)
    return order


def score_largest_gap(edges, ctx, rng):
    """Rank by the unconditional slack that conditioning could remove."""
    scored = []
    for e in edges:
        k = ctx["index"][e]
        gap = 0.0
        for node in e:
            for other_k in range(ctx["d_hid"]):
                other = e[1] if node == e[0] else e[0]
                gap = max(gap, ctx["z1_hi"][other][other_k] - ctx["z1_lo"][other][other_k])
        scored.append((gap, e))
    scored.sort(key=lambda t: (-t[0], t[1]))
    return [e for _, e in scored]


def score_measured_gain(edges, ctx, rng):
    """TEACHER. Actually runs the conditional pass for every edge and ranks by the
    measured tightening. Exact but costs one conditional propagation per edge; the cost
    is counted and reported, never hidden."""
    gains = []
    for e in edges:
        cond = ctx["conditional_z1"](e, count_as="teacher")
        if cond is None:
            gains.append((float("inf"), e))          # provably impossible: free win
            continue
        g = 0.0
        for u in range(ctx["n"]):
            for k in range(ctx["d_hid"]):
                g += max(0.0, ctx["z1_hi"][u][k] - cond["hi"][u][k])
        gains.append((g, e))
    gains.sort(key=lambda t: (-t[0], t[1]))
    return [e for _, e in gains]


SCORERS = {"random": score_random, "largest_gap": score_largest_gap,
           "measured_gain": score_measured_gain}


# ------------------------------------------------------------------------- the bounder
class ConditionalBounder(ResidualBounder):
    def __init__(self, model, X, base, oracle, forced=(), forbidden=(),
                 budget: int = 8, scorer: str = "largest_gap", seed: int = 0,
                 conc_mode: str = "family", cache: BoundCache | None = None,
                 alpha_fn=None, policy=None):
        super().__init__(model, X, base)
        self.oracle = oracle
        self.forced = tuple(sorted(tuple(e) for e in forced))
        self.forbidden = tuple(sorted(tuple(e) for e in forbidden))
        self.budget = int(budget)
        self.scorer_name = scorer
        self.rng = random.Random(seed)
        self.conc_mode = conc_mode
        self.cache = cache if cache is not None else BoundCache()
        self.alpha_fn = alpha_fn
        self.policy = policy
        self.counter = {"oracle_calls": 0}
        self.model_hash = _model_hash(model)
        self.family_hash = oracle.feasibility(self.forced, self.forbidden).family_hash
        self.index = {e: k for k, e in enumerate(self.cands)}
        self.report = {"conditional_passes": 0, "teacher_passes": 0,
                       "proved_impossible": 0, "conditional_used": 0,
                       "fallback_used": 0, "budget": self.budget,
                       "scorer": scorer}

    # -- concretiser -------------------------------------------------------
    def _make_conc(self, forced, forbidden):
        if self.conc_mode == "box":
            return BOX
        return FamilyConc(self.oracle, self.cands, forced, forbidden, self.counter)

    # -- earlier layers, optionally under a condition ----------------------
    def _upto_z1(self, forced, forbidden, cond=()):
        """Recursion target: M1 -> H1 -> Z1 under the given mask.

        Layer 1 is exactly affine in z, so the only approximation before Z1 is the ReLU.
        Concretising with the conditioned family is what makes the result tighter.
        """
        key = CacheKey.make(self.model_hash, self.family_hash, forced, forbidden,
                            "Z1", -1, -1, cond)
        hit = self.cache.get(key)
        if hit is not None:
            return hit
        conc = self._make_conc(forced, forbidden)
        self.conc = conc
        Z0 = self.X @ self.W0 + self.b0
        M1lo, M1hi = self._message_affine(Z0)
        H1lo, H1hi, kinds = [], [], []
        for i in range(self.n):
            rl, rh, rk = [], [], []
            for k in range(len(M1lo[i])):
                alpha = self.alpha_fn("H1", i, k) if self.alpha_fn else None
                a, b, kind = relu_bounds(M1lo[i][k], M1hi[i][k], conc, alpha)
                rl.append(a); rh.append(b); rk.append(kind)
            H1lo.append(rl); H1hi.append(rh); kinds.append(rk)
        Z1lo, Z1hi = self._linear(H1lo, H1hi, self.W1, self.b1)
        out = {"Z1lo": Z1lo, "Z1hi": Z1hi,
               "lo": [[conc.lo(a) for a in row] for row in Z1lo],
               "hi": [[conc.hi(a) for a in row] for row in Z1hi],
               "relu_kinds": kinds, "conc": conc,
               "M1": (M1lo, M1hi), "H1": (H1lo, H1hi)}
        return self.cache.put(key, out)

    def _conditional_z1(self, e, count_as="budget"):
        """Z1 bounds under `forced + {e}`. None iff e is PROVED impossible.

        Only a genuine cache miss is charged: a repeated query for the same
        (model, family, mask, condition) is free, so the budget limits real recursive
        work rather than call sites.
        """
        forced = tuple(sorted(set(self.forced) | {e}))
        feas = self.oracle.feasibility(forced, self.forbidden)
        self.counter["oracle_calls"] += 1
        if not feas.feasible:
            return None
        key = CacheKey.make(self.model_hash, self.family_hash, forced, self.forbidden,
                            "Z1", -1, -1, (e,))
        if key not in self.cache.store:
            self.report["conditional_passes" if count_as == "budget"
                        else "teacher_passes"] += 1
        return self._upto_z1(forced, self.forbidden, cond=(e,))

    # -- the public conditional support -----------------------------------
    def conditional_message_support(self, node: int, direction: int, edge: tuple,
                                    uncond) -> QBound:
        """Certified q_e bounding the coordinate `direction` of the message carried into
        `node` along `edge`, over valid graphs with edge present."""
        other = edge[1] if node == edge[0] else edge[0]
        u_hi = uncond["hi"][other][direction]
        cond = self._conditional_z1(edge)
        if cond is None:
            self.report["proved_impossible"] += 1
            return QBound(edge, node, direction, float("-inf"), "impossible",
                          u_hi, None, True, self.assumptions())
        c_hi = cond["hi"][other][direction]
        # retain the min of valid bounds: a conditional pass can only help
        if c_hi < u_hi:
            self.report["conditional_used"] += 1
            return QBound(edge, node, direction, c_hi, "conditional", u_hi, c_hi, False,
                          self.assumptions())
        self.report["fallback_used"] += 1
        return QBound(edge, node, direction, u_hi, "unconditional_fallback",
                      u_hi, c_hi, False, self.assumptions())

    def assumptions(self):
        return ("z in [0,1]^E relaxation", "float64, no directed rounding",
                "DIAGNOSTIC until roundoff is bounded",
                f"conc_mode={self.conc_mode}", f"budget={self.budget}")

    # -- edge selection ----------------------------------------------------
    def _select_edges(self, uncond):
        alive = []
        for e in self.cands:
            if e in self.forbidden:
                continue
            alive.append(e)
        ctx = {"index": self.index, "d_hid": len(uncond["hi"][0]), "n": self.n,
               "z1_hi": uncond["hi"], "z1_lo": uncond["lo"],
               "conditional_z1": self._conditional_z1, "teacher_passes": 0}
        if self.policy is not None:
            order = self.policy.rank(alive, ctx, self)
        else:
            order = SCORERS[self.scorer_name](alive, ctx, self.rng)
        return order[:self.budget]

    # -- full propagation --------------------------------------------------
    def propagate(self):
        trace = Trace()
        uncond = self._upto_z1(self.forced, self.forbidden)
        conc = uncond["conc"]
        self.conc = conc
        trace.layers["M1"] = uncond["M1"]
        trace.layers["H1"] = uncond["H1"]
        trace.layers["Z1"] = (uncond["Z1lo"], uncond["Z1hi"])
        trace.relu_kinds["H1"] = uncond["relu_kinds"]

        chosen = set(self._select_edges(uncond))
        qbounds = {}

        d_hid = len(uncond["hi"][0])
        M2lo, M2hi = [], []
        for i in range(self.n):
            row_lo, row_hi = [], []
            for k in range(d_hid):
                acc_lo = Aff.constant(0.0, self.m)
                acc_hi = Aff.constant(0.0, self.m)
                for j in range(self.n):
                    if self.B[i, j] != 0.0:
                        acc_lo = acc_lo + uncond["Z1lo"][j][k] * self.B[i, j]
                        acc_hi = acc_hi + uncond["Z1hi"][j][k] * self.B[i, j]
                for (e_idx, other) in self.incident[i]:
                    e = self.cands[e_idx]
                    lowers, uppers = mccormick(e_idx, uncond["Z1lo"][other][k],
                                               uncond["Z1hi"][other][k], self.m, conc)
                    trace.envelope_terms += 4
                    if e in chosen:
                        q = self.conditional_message_support(i, k, e, uncond)
                        qbounds[(i, k, e)] = q
                        if q.proved_impossible:
                            uppers = uppers + [Aff.constant(0.0, self.m)]
                            lowers = lowers + [Aff.constant(0.0, self.m)]
                        else:
                            uppers = uppers + [Aff.unit(e_idx, self.m, q.value)]
                    acc_lo = acc_lo + _tightest_lower(lowers, conc)
                    acc_hi = acc_hi + _tightest_upper(uppers, conc)
                row_lo.append(acc_lo); row_hi.append(acc_hi)
            M2lo.append(row_lo); M2hi.append(row_hi)
        trace.layers["M2"] = (M2lo, M2hi)

        H2lo, H2hi, kinds2 = [], [], []
        for i in range(self.n):
            rl, rh, rk = [], [], []
            for k in range(d_hid):
                alpha = self.alpha_fn("H2", i, k) if self.alpha_fn else None
                a, b, kind = relu_bounds(M2lo[i][k], M2hi[i][k], conc, alpha)
                rl.append(a); rh.append(b); rk.append(kind)
            H2lo.append(rl); H2hi.append(rh); kinds2.append(rk)
        trace.layers["H2"] = (H2lo, H2hi)
        trace.relu_kinds["H2"] = kinds2

        Slo = [Aff.constant(0.0, self.m) for _ in range(d_hid)]
        Shi = [Aff.constant(0.0, self.m) for _ in range(d_hid)]
        for i in range(self.n):
            for k in range(d_hid):
                Slo[k] = Slo[k] + H2lo[i][k]
                Shi[k] = Shi[k] + H2hi[i][k]
        trace.layers["S"] = ([Slo], [Shi])

        wp, wn = _split(self.wr)
        out_lo = Aff.constant(self.br, self.m)
        out_hi = Aff.constant(self.br, self.m)
        for k in range(d_hid):
            if wp[k] != 0.0:
                out_lo = out_lo + Slo[k] * wp[k]
                out_hi = out_hi + Shi[k] * wp[k]
            if wn[k] != 0.0:
                out_lo = out_lo + Shi[k] * wn[k]
                out_hi = out_hi + Slo[k] * wn[k]
        trace.layers["out"] = ([[out_lo]], [[out_hi]])

        self.report.update(self.counter)
        self.report["cache"] = self.cache.stats()
        self.report["qbounds"] = len(qbounds)
        self.report["chosen_edges"] = sorted(chosen)
        return out_lo, out_hi, trace, qbounds


# ------------------------------------------------------------------- the final bound
def conditional_upper_bound(model, X, oracle, base=None, forced=(), forbidden=(),
                            budget=8, scorer="largest_gap", seed=0,
                            conc_mode="family", cache=None, alpha_fn=None, policy=None):
    """U = b + gamma*d + support(a + gamma*c), with the residual bounded conditionally.

    The affine residual term is folded into the objective BEFORE the support call, so a
    single oracle query maximises the combined coefficient vector over the family.
    """
    from .bounds import ASSUMPTIONS, BoundResult
    from .model import backbone_adjacency
    if base is None:
        base = backbone_adjacency(model.n_nodes, dtype=model.dtype)
    feas = oracle.feasibility(forced, forbidden)
    a_dict = model.additive_coefficients(X)
    b, gamma = model.anchor_bias(X), model.gamma
    m = len(model.candidates)

    if not feas.feasible:
        return BoundResult(float("-inf"), float("-inf"), (), 0.0, "conditional",
                           FLOAT_UNVERIFIED, INFEASIBLE, feas.family_hash,
                           _model_hash(model), ASSUMPTIONS, {"reason": feas.reason})

    report = {}
    if not model.residual_active:
        c_vec, d = np.zeros(m), 0.0
    else:
        cb = ConditionalBounder(model, X, base, oracle, forced, forbidden, budget,
                                scorer, seed, conc_mode, cache, alpha_fn, policy)
        lo_aff, hi_aff, trace, qb = cb.propagate()
        use = hi_aff if gamma >= 0 else lo_aff
        c_vec, d = use.coef, use.const
        report = dict(cb.report)
        for kinds in trace.relu_kinds.values():
            for row in kinds:
                for kd in row:
                    report["relu_" + kd] = report.get("relu_" + kd, 0) + 1

    comb = {e: float(a_dict[e] + gamma * c_vec[k])
            for k, e in enumerate(model.candidates)}
    r = oracle.support(comb, forced, forbidden)
    upper = b + gamma * d + r.value
    witness = r.witness
    report["final_oracle_calls"] = report.get("oracle_calls", 0) + 1
    report["conditional_upper"] = float(upper)

    # RETAIN THE MIN OF VALID BOUNDS. The conditional pipeline concretises over the
    # family while the inherited baseline concretises over the box; the two relaxations
    # are both sound but are NOT pointwise comparable, so the conditional route can lose
    # on an individual instance. Keeping the minimum makes the conditional bound
    # dominate the fallback by construction rather than by hope.
    from .bounds import certified_upper_bound as _baseline
    fb = _baseline(model, X, oracle, base, forced, forbidden, mode="combined")
    report["fallback_upper"] = float(fb.upper)
    if fb.upper < upper:
        report["used"] = "unconditional_fallback"
        upper, witness = fb.upper, fb.incumbent_witness
    else:
        report["used"] = "conditional"

    with torch.no_grad():
        inc = float(model(X, model.indicator(witness), base))
    return BoundResult(float(upper), inc, tuple(witness), float(upper - inc),
                       "conditional", FLOAT_UNVERIFIED, PROVED, feas.family_hash,
                       _model_hash(model), ASSUMPTIONS, report)
