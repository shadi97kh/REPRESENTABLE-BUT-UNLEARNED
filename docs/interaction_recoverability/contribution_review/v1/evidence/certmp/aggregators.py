"""Aggregator registry with the property that actually matters.

The paper's original table listed four aggregators and two verdicts. The real condition
is more general: endpoint exactness holds iff the aggregator is MONOTONE UNDER MULTISET
INCLUSION, in either direction.

  isotone   A(S) <= A(T) whenever S subset-of T   -> min at empty set, max at full set
  antitone  A(S) >= A(T) whenever S subset-of T   -> max at empty set, min at full set
  neither                                          -> interior optimum possible, 2 passes unsound

Both monotone directions give exact endpoints, so both are certifiable; only the roles
of the two endpoints swap. That is why `min` belongs in the certified class even though
adding edges lowers its output, and why `logsumexp` belongs even though it appears
nowhere in the original four-row table.
"""
import numpy as np

NEG_INF, POS_INF = -1e9, 1e9

def agg_sum(A, Z):     return A @ Z
def agg_mean(A, Z):    return (A @ Z) / A.sum(1, keepdims=True).clip(min=1)
def agg_degnorm(A, Z):
    d = A.sum(1).clip(min=1) ** -0.5
    return (d[:, None] * A * d[None, :]) @ Z
def agg_max(A, Z):     return np.where(A[:, :, None] > 0, Z[None, :, :], NEG_INF).max(axis=1)
def agg_min(A, Z):     return np.where(A[:, :, None] > 0, Z[None, :, :], POS_INF).min(axis=1)
def agg_logsumexp(A, Z):
    M = np.where(A[:, :, None] > 0, Z[None, :, :], NEG_INF)
    mx = M.max(axis=1, keepdims=True)
    return (mx + np.log(np.exp(M - mx).sum(axis=1, keepdims=True)))[:, 0, :]
def agg_std(A, Z):
    W = (A > 0).astype(float)[:, :, None]
    cnt = W.sum(1).clip(min=1)
    mu = (W * Z[None, :, :]).sum(1) / cnt
    var = (W * (Z[None, :, :] - mu[:, None, :]) ** 2).sum(1) / cnt
    return np.sqrt(var + 1e-12)

# name -> (fn, predicted monotone direction under multiset inclusion)
AGGREGATORS = {
    "sum":        (agg_sum,        "isotone"),
    "max":        (agg_max,        "isotone"),
    "logsumexp":  (agg_logsumexp,  "isotone"),
    "min":        (agg_min,        "antitone"),
    "mean":       (agg_mean,       "neither"),
    "degnorm":    (agg_degnorm,    "neither"),
    "std":        (agg_std,        "neither"),
}
CERTIFIABLE = tuple(n for n, (_, d) in AGGREGATORS.items() if d in ("isotone", "antitone"))

def direction(name):  return AGGREGATORS[name][1]
def apply(name, A, Z): return AGGREGATORS[name][0](A, Z)
