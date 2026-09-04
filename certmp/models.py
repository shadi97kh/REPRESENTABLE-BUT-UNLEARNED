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
        # Affine readout. A NON-NEGATIVE scale is a monotone non-decreasing map, and a
        # constant shift is order-preserving, so neither disturbs the extremality theorem.
        # They exist so a trained model can reach the target scale; H5 still needs scale>=0.
        self.scale, self.shift = 1.0, 0.0

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
        return float(self.scale * (H.sum(0) @ o) + self.shift)

    def structurally_certifiable(self):
        return self.nonneg and self.agg in STRUCTURALLY_SAFE and self.scale >= 0


# ---------------------------------------------------------------------------
# Torch trainable version. The certificate is always computed by the numpy
# MonoMPNN above, so the two forwards MUST agree bit-for-bit-ish or the
# certificate does not describe the trained model. export_numpy() + the parity
# test in tests/test_sanity.py enforce that.
#
# Parameterisation mirrors MonoMPNN exactly:
#   numpy stores W raw and applies softplus  -> torch exports W_raw unchanged
#   numpy uses B directly and needs B >= 0   -> torch exports softplus(b_raw)
#   numpy stores out raw and applies softplus-> torch exports out_raw unchanged
# ---------------------------------------------------------------------------
class TorchMonoMPNN:
    """Lazily imports torch so the reachability core stays numpy-only."""

    def __new__(cls, *a, **kw):
        import torch, torch.nn as nn

        class _Impl(nn.Module):
            def __init__(self, d_in, d_hid=8, n_layers=2, agg="sum", seed=0, nonneg=True):
                super().__init__()
                if agg not in STRUCTURALLY_SAFE:
                    raise ValueError(f"agg={agg} is outside the certifiable class {STRUCTURALLY_SAFE}")
                # nonneg=False is the UNCERTIFIABLE matched baseline: identical in every
                # other respect, so the accuracy difference isolates the sign constraint.
                self.nonneg = nonneg
                g = torch.Generator().manual_seed(seed)
                dims = [d_in] + [d_hid] * n_layers
                self.W = nn.ParameterList([
                    nn.Parameter(torch.randn(dims[i], dims[i + 1], generator=g) * 0.6)
                    for i in range(n_layers)])
                self.b = nn.ParameterList([
                    nn.Parameter(torch.randn(dims[i + 1], generator=g).abs() * 0.1)
                    for i in range(n_layers)])
                self.out = nn.Parameter(torch.randn(dims[-1], generator=g) * 0.6)
                # softplus(log_scale) >= 0 keeps the readout monotone; shift is free
                self.log_scale = nn.Parameter(torch.tensor(0.5413))   # softplus -> ~1.0
                self.shift = nn.Parameter(torch.zeros(()))
                self.agg, self.d_in, self.d_hid, self.n_layers = agg, d_in, d_hid, n_layers

            def set_output_affine(self, scale, shift):
                """Put the readout on the target's scale before training starts."""
                import math
                with torch.no_grad():
                    s = max(float(scale), 1e-12)
                    self.log_scale.fill_(math.log(math.expm1(s)) if s < 20 else s)
                    self.shift.fill_(float(shift))

            def forward(self, X, A):
                """X (B,N,F) non-negative, A (B,N,N) with self-loops. Returns (B,)."""
                H = X
                for li, W in enumerate(self.W):
                    sp_ = torch.nn.functional.softplus
                    Wl = sp_(W) if self.nonneg else W
                    bl = sp_(self.b[li]) if self.nonneg else self.b[li]
                    Z = H @ Wl + bl
                    if self.agg == "sum":
                        M = torch.bmm(A, Z)
                    else:  # max
                        M = torch.where(A.unsqueeze(-1) > 0, Z.unsqueeze(1),
                                        torch.full_like(Z.unsqueeze(1), -1e9)).max(dim=2).values
                    H = torch.relu(M)
                sp = torch.nn.functional.softplus
                o = sp(self.out) if self.nonneg else self.out
                return sp(self.log_scale) * (H.sum(1) @ o) + self.shift

            def export_numpy(self):
                """Hand the trained weights to the certifiable numpy model."""
                import numpy as np
                m = MonoMPNN(self.d_in, self.d_hid, self.n_layers, agg=self.agg,
                             nonneg=self.nonneg, seed=0)
                sp = torch.nn.functional.softplus
                m.W = [w.detach().cpu().numpy().astype(float) for w in self.W]
                m.B = [(sp(b) if self.nonneg else b).detach().cpu().numpy().astype(float)
                       for b in self.b]
                m.out = self.out.detach().cpu().numpy().astype(float)
                m.scale = float(sp(self.log_scale).detach().cpu())
                m.shift = float(self.shift.detach().cpu())
                if self.nonneg:
                    assert all(b.min() >= 0 for b in m.B), "H2 violated: negative bias after export"
                    assert m.scale >= 0, "H5 violated: negative output scale after export"
                return m

        return _Impl(*a, **kw)
