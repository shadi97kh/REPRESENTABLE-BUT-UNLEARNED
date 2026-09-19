"""Ensemble readouts: yhat = b + a^T mu + gamma * E[g(X, z)].

WHY THIS FORM. The measured label is a readout of the whole ensemble, not of one latent
structure. Taking the expectation of the model over the prior,

    E_p[f(X, z)] = b(X) + a(X)^T E_p[z] + gamma * E_p[g(X, z)]
                 = b(X) + a(X)^T mu     + gamma * E_p[g(X, z)]

by linearity of the additive part. So the anchor needs only FIRST MOMENTS mu, and only
the nonlinear residual needs an expectation. No latent structure is ever labelled with
the ensemble outcome.

TWO PRIORS, WITH DIFFERENT ENTITLEMENTS.

  `declared`  a Gibbs measure over the validated non-crossing grammar with declared
              per-pair coefficients. Its marginals are EXACT (the oracle's partition
              function is checked against exhaustive enumeration), and sampling is exact
              sequential conditioning, so every sample is a valid member by construction.

  `vienna`    ViennaRNA base-pair probabilities. A different grammar, not validated here,
              so its marginals are reported as `unvalidated_grammar` and no exactness is
              claimed for them.

Exact marginals are used only for the validated grammar. Nonlinear expectations are
always sampled, with the standard error recorded alongside.

DETERMINISTIC MEAN BRACKET. The affine envelopes satisfy lo(z) <= g(z) <= hi(z)
pointwise, so taking expectations gives

    lo.coef^T mu + lo.const  <=  E[g]  <=  hi.coef^T mu + hi.const

which needs only mu -- no sampling, no Monte Carlo error. It brackets the same quantity
the sampler estimates, so the two are cross-checked against each other.

CHEMISTRY. Only the native target window is ever folded. Guide chemistry lives in X.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
import torch

from ..bounds import ResidualBounder

VALIDATED = "validated_grammar_exact"
UNVALIDATED = "unvalidated_grammar"


@dataclass
class Marginals:
    mu: np.ndarray
    status: str
    grammar: str
    normalisation_ok: bool
    notes: dict = field(default_factory=dict)


@dataclass
class Expectation:
    mean: float
    stderr: float
    n_samples: int
    n_valid: int
    method: str
    bracket: tuple | None = None
    within_bracket: bool | None = None


# ---------------------------------------------------------------- declared prior
class DeclaredGibbsPrior:
    """exp(<theta, z>) / Z over the oracle's validated family. Fixed and declared."""

    def __init__(self, oracle, theta: dict, name="declared"):
        self.oracle = oracle
        self.theta = {tuple(k): float(v) for k, v in theta.items()}
        self.name = name
        self.grammar = f"{oracle.family_id}:non_crossing_canonical"

    def coefficients(self):
        return {e: self.theta.get(e, 0.0) for e in self.oracle.ground_set()}

    def marginals(self, forced=(), forbidden=()) -> Marginals:
        r = self.oracle.marginals(self.coefficients(), forced, forbidden)
        g = self.oracle.ground_set()
        mu = np.array([r.value.get(e, 0.0) for e in g]) if r.value else np.zeros(len(g))
        ok = bool(np.all(mu >= -1e-12) and np.all(mu <= 1 + 1e-12))
        return Marginals(mu, VALIDATED, self.grammar, ok,
                         {"partition_calls": r.stats.get("partition_calls")})

    def sample(self, n, rng, forced=(), forbidden=(), method="backtrace"):
        """Draw exactly from the Gibbs measure over the family.

        `backtrace`  one O(n^3) inside pass, then each draw is a single top-down walk
                     that picks a production with probability proportional to its inside
                     weight. This is the standard stochastic backtrace and is exact.
        `sequential` decide one edge at a time from partition-function ratios. Also
                     exact, needs no chart access, and costs |E| partition evaluations
                     per draw -- kept as an independent cross-check of `backtrace`, not
                     as the production path.

        Either way every draw is a valid member: `backtrace` only ever takes productions
        the recurrence itself admits, and `sequential` only ever conditions on feasible
        masks.
        """
        if method == "sequential" or not hasattr(self.oracle, "log_inside"):
            return self._sample_sequential(n, rng, forced, forbidden)
        return self._sample_backtrace(n, rng, forced, forbidden)

    def _sample_backtrace(self, n, rng, forced=(), forbidden=()):
        o = self.oracle
        coef = self.coefficients()
        chart = o.log_inside(coef, forced, forbidden)
        W = chart["W"]
        if W is None or chart["value"] == float("-inf"):
            raise ValueError("cannot sample from an empty family")
        forced_s, forbidden_s = set(chart["forced"]), set(chart["forbidden"])
        partner = {}
        for (i, j) in forced_s:
            partner[i] = j; partner[j] = i
        N = o.n

        def walk(i, j, acc):
            while i <= j:
                p = partner.get(i)
                options, weights = [], []
                if p is None:
                    w = W[i + 1][j + 1]
                    if w != float("-inf"):
                        options.append(("unpaired", None)); weights.append(w)
                for k in o._pairable[i]:
                    if k > j:
                        break
                    if p is not None and k != p:
                        continue
                    if (i, k) in forbidden_s:
                        continue
                    q = partner.get(k)
                    if q is not None and q != i:
                        continue
                    inner, outer = W[i + 1][k], W[k + 1][j + 1]
                    if inner == float("-inf") or outer == float("-inf"):
                        continue
                    weights.append(coef.get((i, k), 0.0) + inner + outer)
                    options.append(("pair", k))
                if not options:
                    return
                mx = max(weights)
                probs = [math.exp(w - mx) for w in weights]
                total = sum(probs)
                r = rng.random() * total
                pick = 0
                for idx, pr in enumerate(probs):
                    r -= pr
                    if r <= 0:
                        pick = idx; break
                kind, k = options[pick]
                if kind == "unpaired":
                    i += 1
                    continue
                acc.append((i, k))
                walk(i + 1, k - 1, acc)
                i = k + 1

        out = []
        for _ in range(n):
            acc = []
            walk(0, N - 1, acc)
            out.append(tuple(sorted(acc)))
        return out

    def _sample_sequential(self, n, rng, forced=(), forbidden=()):
        coef = self.coefficients()
        ground = list(self.oracle.ground_set())
        out = []
        for _ in range(n):
            f, b = set(forced), set(forbidden)
            for e in ground:
                if e in f or e in b:
                    continue
                z_all = self.oracle.log_partition(coef, tuple(sorted(f)),
                                                  tuple(sorted(b)))
                z_on = self.oracle.log_partition(coef, tuple(sorted(f | {e})),
                                                 tuple(sorted(b)))
                if z_on.value == float("-inf"):
                    b.add(e)
                    continue
                p = math.exp(z_on.value - z_all.value)
                (f if rng.random() < p else b).add(e)
            out.append(tuple(sorted(f)))
        return out


# ------------------------------------------------------------------ vienna prior
class ViennaPrior:
    """ViennaRNA base-pair probabilities. A DIFFERENT grammar, not validated here."""

    def __init__(self, sequence, oracle, name="vienna"):
        self.sequence = sequence
        self.oracle = oracle
        self.name = name
        self.grammar = "ViennaRNA:turner_nearest_neighbour"

    def marginals(self, forced=(), forbidden=()) -> Marginals:
        import RNA
        fc = RNA.fold_compound(self.sequence)
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp = fc.bpp()
        g = self.oracle.ground_set()
        mu = np.array([bpp[i + 1][j + 1] for (i, j) in g])
        ok = bool(np.all(mu >= -1e-12) and np.all(mu <= 1 + 1e-12))
        return Marginals(mu, UNVALIDATED, self.grammar, ok,
                         {"note": "grammar not validated in this repository; no "
                                  "exactness is claimed for these marginals"})

    def sample(self, n, rng, forced=(), forbidden=()):
        import RNA
        RNA.cvar.uniq_ML = 1
        RNA.init_rand(rng.randrange(1 << 30))
        fc = RNA.fold_compound(self.sequence)
        _, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        allowed = set(self.oracle.ground_set())
        out = []
        for s in fc.pbacktrack(n):
            stack, pairs = [], []
            for i, c in enumerate(s):
                if c == "(":
                    stack.append(i)
                elif c == ")" and stack:
                    pairs.append((stack.pop(), i))
            keep = tuple(sorted(p for p in pairs if p in allowed))
            out.append(keep)
        return out


# ------------------------------------------------------------------- expectations
def expected_residual(model, X, base, prior, n_samples, rng,
                      forced=(), forbidden=()) -> Expectation:
    """Monte Carlo E[g], with the standard error recorded, never dropped."""
    if not model.residual_active:
        return Expectation(0.0, 0.0, 0, 0, "residual_inactive")
    samples = prior.sample(n_samples, rng, forced, forbidden)
    vals, valid = [], 0
    for s in samples:
        if not prior.oracle.feasibility(forced=s).feasible:
            continue
        valid += 1
        z = model.indicator(s)
        with torch.no_grad():
            vals.append(float(model.g(X, model.adjacency(z, base))))
    if not vals:
        return Expectation(float("nan"), float("nan"), len(samples), 0, "no_valid_sample")
    arr = np.array(vals)
    se = float(arr.std(ddof=1) / math.sqrt(len(arr))) if len(arr) > 1 else float("nan")
    return Expectation(float(arr.mean()), se, len(samples), valid, "monte_carlo")


def mean_bracket(model, X, base, mu: np.ndarray) -> tuple:
    """Deterministic bracket on E[g] from the affine envelopes and mu alone.

    Sound because lo(z) <= g(z) <= hi(z) holds pointwise on the family, and expectation
    is monotone. Requires only first moments; no sampling error enters.
    """
    if not model.residual_active:
        return (0.0, 0.0)
    lo, hi, _ = ResidualBounder(model, X, base).propagate()
    return (float(lo.coef @ mu + lo.const), float(hi.coef @ mu + hi.const))


def predict_ensemble(model, X, mu, e_g) -> float:
    """yhat = b + a^T mu + gamma * E[g]."""
    b = model.anchor_bias(X)
    a = model.additive_coefficients(X)
    avec = np.array([a[e] for e in model.candidates])
    return float(b + avec @ mu + model.gamma * e_g)


def ensemble_bracket(model, X, mu, bracket) -> tuple:
    """Deterministic bracket on yhat, from the bracket on E[g]."""
    b = model.anchor_bias(X)
    a = model.additive_coefficients(X)
    avec = np.array([a[e] for e in model.candidates])
    lin = float(b + avec @ mu)
    g = model.gamma
    lo, hi = (bracket if g >= 0 else (bracket[1], bracket[0]))
    return (lin + g * lo, lin + g * hi)
