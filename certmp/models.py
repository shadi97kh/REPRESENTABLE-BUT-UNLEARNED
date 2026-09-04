"""Monotone MPNN. numpy reference model (no torch needed) + torch trainable version.

Monotonicity hypotheses, ALL required for the structural theorem:
  H1 non-negative message/update weights   (softplus reparameterisation)
  H2 non-negative biases
  H3 monotone non-decreasing activation
  H4 non-negative node features X >= 0
  H5 non-negative readout
"""
import numpy as np

def softplus(z): return np.log1p(np.exp(-np.abs(z))) + np.maximum(z, 0)
RELU = lambda z: np.maximum(z, 0)
AGGS = ("sum", "mean", "degnorm", "max")
STRUCTURALLY_SAFE = ("sum", "max")     # the theorem's scope

class MonoMPNN:
    def __init__(self, d_in, d_hid=8, n_layers=2, agg="sum", nonneg=True, seed=0):
        r = np.random.RandomState(seed); dims = [d_in] + [d_hid]*n_layers
        self.W = [r.randn(dims[i], dims[i+1])*0.6 for i in range(n_layers)]
        self.B = [np.abs(r.randn(dims[i+1]))*0.1 for i in range(n_layers)]
        self.out = r.randn(dims[-1])*0.6
        self.agg, self.nonneg = agg, nonneg

    def _w(self, M): return softplus(M) if self.nonneg else M

    def forward(self, X, A):
        H = np.asarray(X, dtype=float)
        for li, W in enumerate(self.W):
            Z = H @ self._w(W) + self.B[li]
            if self.agg == "sum":
                M = A @ Z
            elif self.agg == "mean":
                M = (A @ Z) / A.sum(1, keepdims=True).clip(min=1)
            elif self.agg == "degnorm":
                d = A.sum(1).clip(min=1) ** -0.5
                M = (d[:, None] * A * d[None, :]) @ Z
            elif self.agg == "max":
                M = np.where(A[:, :, None] > 0, Z[None, :, :], -1e9).max(axis=1)
            else:
                raise ValueError(self.agg)
            H = RELU(M)
        o = softplus(self.out) if self.nonneg else self.out
        return float(H.sum(0) @ o)

    def structurally_certifiable(self):
        return self.nonneg and self.agg in STRUCTURALLY_SAFE
