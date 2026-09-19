"""Exact reachability. Two forward passes vs 2^k brute force."""
import itertools, numpy as np
from .void import check_live, check_nonempty_lattice, VoidRun


def build_A(n, mandatory, optional, keep, self_loops=True):
    A = np.eye(n) if self_loops else np.zeros((n, n))
    for (i, j) in mandatory: A[i, j] = A[j, i] = 1.0
    for k, (i, j) in enumerate(optional):
        if keep[k]: A[i, j] = A[j, i] = 1.0
    return A


def _endpoints(model, n, mandatory, optional, X):
    """Both lattice endpoints, ordered (lo, hi) for this model's monotone direction.
    Isotone: min at the empty optional set, max at the full one. Antitone: reversed."""
    y_empty = model.forward(X, build_A(n, mandatory, optional, [0] * len(optional)))
    y_full = model.forward(X, build_A(n, mandatory, optional, [1] * len(optional)))
    antitone = getattr(model, "endpoint_direction", lambda: "isotone")() == "antitone"
    return (y_full, y_empty) if antitone else (y_empty, y_full)


def exact_interval(model, n, mandatory, optional, X):
    """THE THEOREM. Valid only when model.structurally_certifiable().

    The reachable extremes are attained at the two lattice endpoints. Which endpoint
    carries which extreme depends on the aggregator's monotone direction, so both isotone
    and antitone aggregators are certifiable in two passes.
    """
    if not model.structurally_certifiable():
        from .aggregators import CERTIFIABLE
        from .models import MONOTONE_ACTS
        raise ValueError(
            f"agg={model.agg} act={getattr(model, 'act', 'relu')} nonneg={model.nonneg}: "
            f"outside the certified class. Need non-negative parameters, an aggregator "
            f"monotone under multiset inclusion {CERTIFIABLE}, and a monotone activation "
            f"{MONOTONE_ACTS}.")
    if np.asarray(X).min() < 0:
        raise ValueError("hypothesis H4 violated: node features must be non-negative")
    lo, hi = _endpoints(model, n, mandatory, optional, X)
    return lo, hi, 2                      # 2 = forward passes used


def brute_force_interval(model, n, mandatory, optional, X, max_k=20):
    """Ground-truth oracle. Exponential; used only to VALIDATE the theorem."""
    k = len(optional)
    check_nonempty_lattice(optional)
    if k > max_k: raise ValueError(f"k={k} too large to enumerate (2^{k})")
    vals = [model.forward(X, build_A(n, mandatory, optional, bits))
            for bits in itertools.product([0, 1], repeat=k)]
    check_live(vals)
    return min(vals), max(vals), 2**k


def endpoint_gap(model, n, mandatory, optional, X):
    """Relative amount by which the TRUE optimum beats the lattice endpoint.
    0.0 => the theorem holds on this instance."""
    blo, bhi, n_states = brute_force_interval(model, n, mandatory, optional, X)
    elo, ehi = _endpoints(model, n, mandatory, optional, X)
    s = max(abs(bhi), abs(blo), 1e-12)
    return {"gap_hi": (bhi - ehi) / s, "gap_lo": (elo - blo) / s,
            "exact": abs(bhi - ehi) / s <= 1e-12 and abs(elo - blo) / s <= 1e-12,
            "n_states": n_states, "true_lo": blo, "true_hi": bhi}
