"""Learned components, each separately switchable.

HARD CONSTRAINT: no network here predicts a bound. Two things are learned, and both are
choices among options that are already valid whatever the network says:

  LearnedSlopePolicy   emits the ReLU lower-bound slope alpha. Lemma 3 in
                       docs/bounds_proof.md holds for EVERY alpha in [0,1], and the
                       value is squashed by a sigmoid and clamped again in
                       `relu_bounds`. A badly trained policy loses tightness; it cannot
                       lose soundness.

  LearnedQueryScorer   emits a RANKING over candidate edges to spend the conditional
                       budget on. The bound for any subset is sound, so the ranking
                       changes which sound bound is computed, never whether it is sound.

Neither returns a number that is used as an upper bound, and the audit re-checks every
bound produced with a policy attached against exhaustive enumeration.
"""
from __future__ import annotations

import numpy as np
import torch
import torch.nn as nn


def edge_features(edges, ctx, n_nodes) -> torch.Tensor:
    """Cheap, model-agnostic descriptors of a candidate edge.

    Deliberately excludes anything that would let the scorer reconstruct a bound: it
    sees slack widths and geometry, not the bound itself.
    """
    d_hid = ctx["d_hid"]
    rows = []
    for (i, j) in edges:
        span = (j - i) / max(n_nodes, 1)
        gap_i = max(ctx["z1_hi"][i][k] - ctx["z1_lo"][i][k] for k in range(d_hid))
        gap_j = max(ctx["z1_hi"][j][k] - ctx["z1_lo"][j][k] for k in range(d_hid))
        mean_i = float(np.mean([ctx["z1_hi"][i][k] - ctx["z1_lo"][i][k]
                                for k in range(d_hid)]))
        mean_j = float(np.mean([ctx["z1_hi"][j][k] - ctx["z1_lo"][j][k]
                                for k in range(d_hid)]))
        centrality = (min(i, n_nodes - 1 - j)) / max(n_nodes, 1)
        rows.append([span, gap_i, gap_j, mean_i, mean_j, centrality,
                     i / max(n_nodes, 1), j / max(n_nodes, 1)])
    return torch.tensor(rows, dtype=torch.float64)


N_EDGE_FEATURES = 8


class LearnedQueryScorer(nn.Module):
    """Ranks edges for the conditional budget. Output is an ordering, never a bound."""

    def __init__(self, d_hid: int = 24, dtype=torch.float64, seed: int | None = None):
        super().__init__()
        if seed is not None:
            torch.manual_seed(seed)
        self.net = nn.Sequential(
            nn.Linear(N_EDGE_FEATURES, d_hid, dtype=dtype), nn.ReLU(),
            nn.Linear(d_hid, d_hid, dtype=dtype), nn.ReLU(),
            nn.Linear(d_hid, 1, dtype=dtype))

    def scores(self, edges, ctx, n_nodes) -> torch.Tensor:
        if not edges:
            return torch.zeros(0, dtype=torch.float64)
        return self.net(edge_features(edges, ctx, n_nodes)).squeeze(-1)

    def rank(self, edges, ctx, bounder):
        with torch.no_grad():
            s = self.scores(list(edges), ctx, bounder.n)
        order = sorted(range(len(edges)), key=lambda t: (-float(s[t]), edges[t]))
        return [edges[t] for t in order]


class LearnedSlopePolicy(nn.Module):
    """Emits alpha in [0,1] per unstable ReLU. Valid for any output by Lemma 3."""

    def __init__(self, d_hid: int = 16, dtype=torch.float64, seed: int | None = None):
        super().__init__()
        if seed is not None:
            torch.manual_seed(seed)
        self.net = nn.Sequential(
            nn.Linear(4, d_hid, dtype=dtype), nn.ReLU(),
            nn.Linear(d_hid, 1, dtype=dtype))
        self._ctx = {}

    def bind(self, layer_stats: dict):
        """layer_stats maps (layer, node, unit) -> (l, u) from the current propagation."""
        self._ctx = layer_stats

    def alpha_fn(self, layer, node, unit):
        lu = self._ctx.get((layer, node, unit))
        if lu is None:
            return None                      # fall back to the deterministic heuristic
        l, u = lu
        if not (l < 0 < u):
            return None
        feats = torch.tensor([[l, u, u - l, u / (u - l)]], dtype=torch.float64)
        with torch.no_grad():
            return float(torch.sigmoid(self.net(feats)).item())


def alpha_search_teacher(evaluate, n_draws: int, rng) -> tuple:
    """Best-of-K random alpha vectors. The teacher for the slope policy.

    `evaluate(alpha_vector) -> upper bound`. Returns (best_alpha, best_upper, n_calls),
    and the call count is the teacher-generation cost for slopes.
    """
    best, best_a = None, None
    for _ in range(n_draws):
        a = rng.random()
        u = evaluate(a)
        if best is None or u < best:
            best, best_a = u, a
    return best_a, best, n_draws
