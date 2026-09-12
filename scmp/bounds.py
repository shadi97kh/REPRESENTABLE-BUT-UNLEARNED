"""Affine bound propagation for f(X, z) = b + a^T z + gamma * g(X, B + P z).

Every intermediate quantity is bracketed by two AFFINE FUNCTIONS OF z,

    lo.coef^T z + lo.const  <=  quantity(z)  <=  hi.coef^T z + hi.const     for all z in [0,1]^E,

so the whole network collapses to one affine upper bound on g, which the support oracle
then maximises over the feasible family in a single call.

WHERE THE BILINEARITY IS. With Z0 = X W0 + b0 constant in z,

    M1 = (B + Pz) Z0        is AFFINE in z          (Z0 does not depend on z)
    M2 = (B + Pz) Z1(z)     is BILINEAR in z        (Z1 depends on z through H1)

so layer 1 needs no relaxation at all and layer 2 needs the product envelope. The same
z_e variable appears in both layers and is NOT duplicated: a fresh copy per layer would
be sound but strictly looser, because it would allow an edge to be present at one depth
and absent at another.

PRODUCT ENVELOPE. For y = e h with e in [0,1] and h in [hl, hu], all four McCormick
inequalities are formed and all four are used:

    y >= e hl                       y >= h + e hu - hu
    y <= h + e hl - hl              y <= e hu

The two containing h are closed by substituting h's affine lower bound where h carries a
positive coefficient in a lower bound, and its affine upper bound in an upper bound. Each
of the four is individually valid, so selecting the tightest per term is sound.

NUMERICAL STATUS. Propagation runs in ordinary float64 with no directed rounding and no
interval arithmetic, so every bound produced here is labelled `float_unverified` and is a
DIAGNOSTIC. It is not a proof. See docs/bounds_proof.md for the real-arithmetic argument
that is being approximated and the plan for making it numerically sound.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import torch

FLOAT_UNVERIFIED = "float_unverified"
EXACT_INTEGER = "exact_integer"
PROVED = "proved"
INFEASIBLE = "infeasible"


# --------------------------------------------------------------------------- affine form
@dataclass
class Aff:
    """coef^T z + const, with z indexed by the candidate universe."""
    coef: np.ndarray
    const: float

    @staticmethod
    def constant(value: float, m: int) -> "Aff":
        return Aff(np.zeros(m), float(value))

    @staticmethod
    def unit(idx: int, m: int, weight: float = 1.0) -> "Aff":
        c = np.zeros(m); c[idx] = weight
        return Aff(c, 0.0)

    def __add__(self, other):
        if isinstance(other, Aff):
            return Aff(self.coef + other.coef, self.const + other.const)
        return Aff(self.coef.copy(), self.const + float(other))

    def __mul__(self, s: float):
        s = float(s)
        return Aff(self.coef * s, self.const * s)

    __rmul__ = __mul__

    def box_min(self) -> float:
        return self.const + float(np.minimum(self.coef, 0.0).sum())

    def box_max(self) -> float:
        return self.const + float(np.maximum(self.coef, 0.0).sum())

    def eval(self, z: np.ndarray) -> float:
        return float(self.coef @ z + self.const)


def _split(W: np.ndarray):
    return np.maximum(W, 0.0), np.minimum(W, 0.0)


class BoxConc:
    """Concretise an affine form over the box [0,1]^E. The default everywhere."""
    def lo(self, a: "Aff") -> float: return a.box_min()
    def hi(self, a: "Aff") -> float: return a.box_max()


BOX = BoxConc()


# ------------------------------------------------------------------------ relaxations
def relu_bounds(lo: Aff, hi: Aff, conc=None, alpha=None):
    """Affine bracket for relu, with the two stable cases handled exactly.

    `conc` concretises the incoming bracket; the default is the box. A tighter
    concretiser (e.g. maximising over the feasible family) only shrinks [l, u], which
    can turn an unstable unit stable and is sound either way.

    `alpha` overrides the lower-bound slope. ANY value in [0, 1] is valid (Lemma 3), so
    a learned alpha cannot make the bound unsound; it is clamped here regardless.
    """
    conc = conc or BOX
    l, u = conc.lo(lo), conc.hi(hi)
    m = len(lo.coef)
    if l >= 0.0:
        return lo, hi, "stable_pos"          # relu is the identity here: exact
    if u <= 0.0:
        return Aff.constant(0.0, m), Aff.constant(0.0, m), "stable_neg"
    # unstable: upper is the chord through (l,0) and (u,u); lower is alpha*x
    s = u / (u - l)
    up = Aff(hi.coef * s, hi.const * s - s * l)
    a = (1.0 if u >= -l else 0.0) if alpha is None else float(np.clip(alpha, 0.0, 1.0))
    low = Aff(lo.coef * a, lo.const * a)
    return low, up, "unstable"


def mccormick(e_idx: int, h_lo: Aff, h_hi: Aff, m: int, conc=None):
    """All four envelope inequalities for y = z_e * h, z_e in [0,1].

    Returns (lower_candidates, upper_candidates), each a list of two affine forms.
    """
    conc = conc or BOX
    hl, hu = conc.lo(h_lo), conc.hi(h_hi)
    lowers = [
        Aff.unit(e_idx, m, hl),                                   # y >= e*hl
        h_lo + Aff.unit(e_idx, m, hu) + (-hu),                    # y >= h + e*hu - hu
    ]
    uppers = [
        h_hi + Aff.unit(e_idx, m, hl) + (-hl),                    # y <= h + e*hl - hl
        Aff.unit(e_idx, m, hu),                                   # y <= e*hu
    ]
    return lowers, uppers


def _tightest_lower(cands, conc=None):
    c = conc or BOX
    return max(cands, key=lambda a: c.lo(a))


def _tightest_upper(cands, conc=None):
    c = conc or BOX
    return min(cands, key=lambda a: c.hi(a))


# ------------------------------------------------------------------------- propagation
@dataclass
class Trace:
    """Every intermediate bracket, so the tests can check each inequality, not only
    the final one."""
    layers: dict = field(default_factory=dict)
    relu_kinds: dict = field(default_factory=dict)
    envelope_terms: int = 0


class ResidualBounder:
    """Affine bounds on g(X, B + P z) as a function of z."""

    def __init__(self, model, X: torch.Tensor, base: torch.Tensor):
        self.model = model
        self.cands = model.candidates
        self.m = len(self.cands)
        self.n = model.n_nodes
        self.X = X.detach().cpu().numpy().astype(float)
        self.B = base.detach().cpu().numpy().astype(float)
        g = model.g
        self.W0 = g.lin0.weight.detach().cpu().numpy().T.astype(float)
        self.b0 = g.lin0.bias.detach().cpu().numpy().astype(float)
        self.W1 = g.lin1.weight.detach().cpu().numpy().T.astype(float)
        self.b1 = g.lin1.bias.detach().cpu().numpy().astype(float)
        self.wr = g.readout.weight.detach().cpu().numpy().reshape(-1).astype(float)
        self.br = float(g.readout.bias.detach().cpu().numpy().reshape(-1)[0])
        self.incident = {i: [] for i in range(self.n)}
        for k, (i, j) in enumerate(self.cands):
            self.incident[i].append((k, j))
            self.incident[j].append((k, i))

    # -- helpers ---------------------------------------------------------
    def _conc(self):
        return getattr(self, "conc", None) or BOX

    def _linear(self, lo, hi, W, bias):
        """Signed linear map applied to an affine bracket, elementwise over nodes."""
        Wp, Wn = _split(W)
        d_out = W.shape[1]
        out_lo, out_hi = [], []
        for i in range(self.n):
            row_lo, row_hi = [], []
            for k in range(d_out):
                acc_lo = Aff.constant(bias[k], self.m)
                acc_hi = Aff.constant(bias[k], self.m)
                for p in range(W.shape[0]):
                    if Wp[p, k] != 0.0:
                        acc_lo = acc_lo + lo[i][p] * Wp[p, k]
                        acc_hi = acc_hi + hi[i][p] * Wp[p, k]
                    if Wn[p, k] != 0.0:
                        acc_lo = acc_lo + hi[i][p] * Wn[p, k]
                        acc_hi = acc_hi + lo[i][p] * Wn[p, k]
                row_lo.append(acc_lo); row_hi.append(acc_hi)
            out_lo.append(row_lo); out_hi.append(row_hi)
        return out_lo, out_hi

    def _message_affine(self, Zc: np.ndarray):
        """M = (B + Pz) Zc with Zc CONSTANT: exactly affine, no relaxation."""
        d = Zc.shape[1]
        lo, hi = [], []
        for i in range(self.n):
            row = []
            for k in range(d):
                a = Aff.constant(float(self.B[i] @ Zc[:, k]), self.m)
                for (e_idx, other) in self.incident[i]:
                    a.coef[e_idx] += Zc[other, k]
                row.append(a)
            lo.append(row); hi.append([Aff(r.coef.copy(), r.const) for r in row])
        return lo, hi

    def _message_bilinear(self, Zlo, Zhi, trace: Trace):
        """M = (B + Pz) Z(z): the B part is a non-negative linear map, the Pz part is
        the product envelope, one term per incident candidate edge."""
        d = len(Zlo[0])
        out_lo, out_hi = [], []
        for i in range(self.n):
            row_lo, row_hi = [], []
            for k in range(d):
                acc_lo = Aff.constant(0.0, self.m)
                acc_hi = Aff.constant(0.0, self.m)
                for j in range(self.n):                      # B >= 0 entrywise
                    if self.B[i, j] != 0.0:
                        acc_lo = acc_lo + Zlo[j][k] * self.B[i, j]
                        acc_hi = acc_hi + Zhi[j][k] * self.B[i, j]
                for (e_idx, other) in self.incident[i]:
                    lowers, uppers = mccormick(e_idx, Zlo[other][k], Zhi[other][k],
                                               self.m, self._conc())
                    trace.envelope_terms += 4
                    acc_lo = acc_lo + _tightest_lower(lowers, self._conc())
                    acc_hi = acc_hi + _tightest_upper(uppers, self._conc())
                row_lo.append(acc_lo); row_hi.append(acc_hi)
            out_lo.append(row_lo); out_hi.append(row_hi)
        return out_lo, out_hi

    def reference_forward(self, z: np.ndarray) -> dict:
        """Exact intermediate values at one z, in numpy from the extracted weights.

        Independent of the affine machinery and of torch: the tests cross-check its
        final value against the torch model, then use its intermediates to check every
        bracket the propagation produced.
        """
        A = self.B.copy()
        for k, (i, j) in enumerate(self.cands):
            A[i, j] += z[k]
            A[j, i] += z[k]
        Z0 = self.X @ self.W0 + self.b0
        M1 = A @ Z0
        H1 = np.maximum(M1, 0.0)
        Z1 = H1 @ self.W1 + self.b1
        M2 = A @ Z1
        H2 = np.maximum(M2, 0.0)
        S = H2.sum(axis=0)
        out = float(self.wr @ S + self.br)
        return {"M1": M1, "H1": H1, "Z1": Z1, "M2": M2, "H2": H2,
                "S": S[None, :], "out": np.array([[out]])}

    # -- the propagation --------------------------------------------------
    def propagate(self) -> tuple[Aff, Aff, Trace]:
        trace = Trace()
        Z0 = self.X @ self.W0 + self.b0                       # constant in z
        M1lo, M1hi = self._message_affine(Z0)
        trace.layers["M1"] = (M1lo, M1hi)

        H1lo, H1hi, kinds = [], [], []
        for i in range(self.n):
            rl, rh, rk = [], [], []
            for k in range(len(M1lo[i])):
                a, b, kind = relu_bounds(M1lo[i][k], M1hi[i][k], self._conc())
                rl.append(a); rh.append(b); rk.append(kind)
            H1lo.append(rl); H1hi.append(rh); kinds.append(rk)
        trace.layers["H1"] = (H1lo, H1hi)
        trace.relu_kinds["H1"] = kinds

        Z1lo, Z1hi = self._linear(H1lo, H1hi, self.W1, self.b1)
        trace.layers["Z1"] = (Z1lo, Z1hi)

        M2lo, M2hi = self._message_bilinear(Z1lo, Z1hi, trace)
        trace.layers["M2"] = (M2lo, M2hi)

        H2lo, H2hi, kinds2 = [], [], []
        for i in range(self.n):
            rl, rh, rk = [], [], []
            for k in range(len(M2lo[i])):
                a, b, kind = relu_bounds(M2lo[i][k], M2hi[i][k], self._conc())
                rl.append(a); rh.append(b); rk.append(kind)
            H2lo.append(rl); H2hi.append(rh); kinds2.append(rk)
        trace.layers["H2"] = (H2lo, H2hi)
        trace.relu_kinds["H2"] = kinds2

        # readout: signed linear on the node-summed hidden state
        d = len(H2lo[0])
        Slo = [Aff.constant(0.0, self.m) for _ in range(d)]
        Shi = [Aff.constant(0.0, self.m) for _ in range(d)]
        for i in range(self.n):
            for k in range(d):
                Slo[k] = Slo[k] + H2lo[i][k]
                Shi[k] = Shi[k] + H2hi[i][k]
        trace.layers["S"] = ([Slo], [Shi])

        wp, wn = _split(self.wr)
        out_lo = Aff.constant(self.br, self.m)
        out_hi = Aff.constant(self.br, self.m)
        for k in range(d):
            if wp[k] != 0.0:
                out_lo = out_lo + Slo[k] * wp[k]
                out_hi = out_hi + Shi[k] * wp[k]
            if wn[k] != 0.0:
                out_lo = out_lo + Shi[k] * wn[k]
                out_hi = out_hi + Slo[k] * wn[k]
        trace.layers["out"] = ([[out_lo]], [[out_hi]])
        return out_lo, out_hi, trace


# ------------------------------------------------------------------------- the bound
@dataclass
class BoundResult:
    upper: float
    incumbent_value: float
    incumbent_witness: tuple
    gap: float
    mode: str
    numerical_status: str
    termination_status: str
    family_hash: str
    model_hash: str
    assumptions: tuple
    stats: dict = field(default_factory=dict)


def _model_hash(model) -> str:
    import hashlib
    h = hashlib.sha256()
    for name, p in sorted(model.state_dict().items()):
        h.update(name.encode())
        h.update(np.ascontiguousarray(p.detach().cpu().numpy()).tobytes())
    h.update(f"gamma={model.gamma}".encode())
    return h.hexdigest()[:16]


ASSUMPTIONS = (
    "z in [0,1]^E relaxation of the family indicator",
    "float64 propagation, no directed rounding, no interval arithmetic",
    "bound is DIAGNOSTIC until roundoff is bounded (docs/bounds_proof.md)",
    "edge variables shared across layers, not duplicated",
)


def certified_upper_bound(model, X, oracle, base=None, forced=(), forbidden=(),
                          mode="combined") -> BoundResult:
    """Upper bound on max_{z in F} f(X, z), with a feasible incumbent.

    mode="combined"      U = b + gamma*d + support(a + gamma*c)   -- one oracle call on
                         the COMBINED coefficient vector.
    mode="terminal_only" U = b + support(a) + gamma*box(g)        -- the separate-maxima
                         baseline, identical predictor weights.

    The incumbent is always a real model evaluation f(X, z*) at a feasible z*, never an
    envelope value.
    """
    from .model import backbone_adjacency
    if base is None:
        base = backbone_adjacency(model.n_nodes, dtype=model.dtype)

    feas = oracle.feasibility(forced, forbidden)
    a_dict = model.additive_coefficients(X)
    b = model.anchor_bias(X)
    gamma = model.gamma
    m = len(model.candidates)
    a_vec = np.array([a_dict[e] for e in model.candidates])

    if not feas.feasible:
        return BoundResult(float("-inf"), float("-inf"), (), 0.0, mode,
                           FLOAT_UNVERIFIED, INFEASIBLE, feas.family_hash,
                           _model_hash(model), ASSUMPTIONS,
                           {"reason": feas.reason})

    stats = {"relu_stable_pos": 0, "relu_stable_neg": 0, "relu_unstable": 0,
             "envelope_inequalities": 0, "oracle_calls": 0}

    if not model.residual_active:
        c_vec = np.zeros(m); d = 0.0
    else:
        lo_aff, hi_aff, trace = ResidualBounder(model, X, base).propagate()
        use = hi_aff if gamma >= 0 else lo_aff        # gamma flips which side is needed
        c_vec, d = use.coef, use.const
        for kinds in trace.relu_kinds.values():
            for row in kinds:
                for k in row:
                    stats["relu_" + k] += 1
        stats["envelope_inequalities"] = trace.envelope_terms

    if mode == "combined":
        comb = {e: float(a_vec[k] + gamma * c_vec[k])
                for k, e in enumerate(model.candidates)}
        r = oracle.support(comb, forced, forbidden)
        stats["oracle_calls"] = 1
        upper = b + gamma * d + r.value
        witness = r.witness
    elif mode == "terminal_only":
        r = oracle.support(a_dict, forced, forbidden)
        stats["oracle_calls"] = 1
        resid_box = (Aff(c_vec, d).box_max() if gamma >= 0
                     else Aff(c_vec, d).box_min())
        upper = b + r.value + gamma * resid_box
        witness = r.witness
    else:
        raise ValueError(f"unknown mode {mode!r}")

    # INCUMBENT: evaluate the true model at a feasible point. Never an envelope value.
    with torch.no_grad():
        inc = float(model(X, model.indicator(witness), base))

    return BoundResult(float(upper), inc, tuple(witness), float(upper - inc), mode,
                       FLOAT_UNVERIFIED, PROVED, feas.family_hash, _model_hash(model),
                       ASSUMPTIONS, stats)
