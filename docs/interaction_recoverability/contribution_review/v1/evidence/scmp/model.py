"""The certified model:  f(X, z) = b(X) + a(X)^T z + gamma * g(X, B + P z)

    b(X), a(X)   the ANCHOR. Depends on X and the candidate universe only.
    z            indicator over the candidate universe; a member of the family.
    B            fixed base adjacency (backbone, self-loops) — never depends on z.
    P            the candidate-to-adjacency map: P z places a symmetric 1 at each
                 selected candidate edge. Same index order as `candidates`.
    g            two-layer signed-weight ReLU MPNN, sum aggregation, signed linear
                 readout. The non-additive part.
    gamma        residual scale. gamma == 0.0 is an EXACT bypass, not a small number.

WHY THE ANCHOR MUST NOT SEE z OR THE MASK. The support oracle optimises a^T z over the
family under many different forced/forbidden masks. If a depended on z, the objective
would change with the point being evaluated and the oracle's optimum would not be an
optimum of f. If a depended on the mask, coefficients would have to be recomputed per
query and the anchor would no longer be the fixed additive object the certificate is
built on. `AnchorEncoder.forward` therefore takes X and the candidate list and nothing
else — enforced by its signature, and tested in tests/model.

DELIBERATELY ABSENT from the first certified model, because bound propagation through
them is not supported here: mean or degree normalisation, attention or any softmax over
neighbours, gating, batch or layer norm. Only sum aggregation appears.
"""
from __future__ import annotations

from typing import Iterable, Sequence

import torch
import torch.nn as nn


def backbone_adjacency(n: int, self_loops: bool = True, chain: bool = True,
                       dtype=torch.float64) -> torch.Tensor:
    """The fixed part B: self-loops and the sequence backbone. Never depends on z."""
    A = torch.eye(n, dtype=dtype) if self_loops else torch.zeros(n, n, dtype=dtype)
    if chain:
        for i in range(n - 1):
            A[i, i + 1] = 1.0
            A[i + 1, i] = 1.0
    return A


class AnchorEncoder(nn.Module):
    """(b, a) from the fixed input. Sees no z, no adjacency, no conditioning mask.

    a_e is a symmetric function of the endpoint attributes, so relabelling the nodes and
    carrying the attributes permutes a consistently and leaves b unchanged.
    """

    def __init__(self, d_in: int, d_hid: int = 16, dtype=torch.float64):
        super().__init__()
        self.d_in = d_in
        self.phi = nn.Sequential(
            nn.Linear(2 * d_in, d_hid, dtype=dtype), nn.ReLU(),
            nn.Linear(d_hid, d_hid, dtype=dtype))
        self.a_out = nn.Linear(d_hid, 1, dtype=dtype)
        self.psi = nn.Sequential(
            nn.Linear(d_in, d_hid, dtype=dtype), nn.ReLU(),
            nn.Linear(d_hid, 1, dtype=dtype))

    def forward(self, X: torch.Tensor, candidates: Sequence[tuple]):
        b = self.psi(X).sum()
        if len(candidates) == 0:
            return b, X.new_zeros(0)
        i_idx = torch.tensor([e[0] for e in candidates], dtype=torch.long)
        j_idx = torch.tensor([e[1] for e in candidates], dtype=torch.long)
        xi, xj = X[i_idx], X[j_idx]
        h = self.phi(torch.cat([xi, xj], dim=-1)) + self.phi(torch.cat([xj, xi], dim=-1))
        return b, self.a_out(h).squeeze(-1)


class SignedResidualMPNN(nn.Module):
    """g: two layers, signed weights, ReLU, sum aggregation, signed linear readout.

    The same adjacency is shared by both layers ("shared edge indicators") — there is no
    per-layer edge parameterisation, so a selected candidate edge means the same thing at
    every depth. No normalisation of any kind is applied to A.
    """

    def __init__(self, d_in: int, d_hid: int = 16, dtype=torch.float64):
        super().__init__()
        self.lin0 = nn.Linear(d_in, d_hid, dtype=dtype)     # signed: unconstrained
        self.lin1 = nn.Linear(d_hid, d_hid, dtype=dtype)
        self.readout = nn.Linear(d_hid, 1, dtype=dtype)     # signed linear readout
        self.act = nn.ReLU()

    def forward(self, X: torch.Tensor, A: torch.Tensor) -> torch.Tensor:
        H = self.act(A @ self.lin0(X))          # layer 1, sum aggregation
        H = self.act(A @ self.lin1(H))          # layer 2, same A
        return self.readout(H.sum(dim=0)).squeeze()


class CertifiedModel(nn.Module):
    """f(X, z) = b(X) + a(X)^T z + gamma * g(X, B + P z)."""

    def __init__(self, d_in: int, candidates: Iterable[tuple], n_nodes: int,
                 d_hid: int = 16, gamma: float = 0.0, use_residual: bool = True,
                 dtype=torch.float64, seed: int | None = None):
        super().__init__()
        if seed is not None:
            torch.manual_seed(seed)
        self.candidates = tuple(tuple(e) for e in candidates)
        self.index = {e: k for k, e in enumerate(self.candidates)}
        self.n_nodes = int(n_nodes)
        self.dtype = dtype
        self.use_residual = bool(use_residual)
        self.gamma = float(gamma)
        self.anchor = AnchorEncoder(d_in, d_hid, dtype=dtype)
        self.g = SignedResidualMPNN(d_in, d_hid, dtype=dtype) if use_residual else None
        # Derived from `candidates`, not learned state, so they are NOT persistent:
        # a state_dict from a model over one candidate ordering must not silently
        # overwrite the index map of a model built over a different ordering.
        ii = [e[0] for e in self.candidates]
        jj = [e[1] for e in self.candidates]
        self.register_buffer("_i", torch.tensor(ii, dtype=torch.long), persistent=False)
        self.register_buffer("_j", torch.tensor(jj, dtype=torch.long), persistent=False)

    # ------------------------------------------------------------------ residual switch
    @property
    def residual_active(self) -> bool:
        """True only when a residual exists AND gamma is not exactly zero."""
        return self.use_residual and self.g is not None and self.gamma != 0.0

    def set_gamma(self, gamma: float):
        self.gamma = float(gamma)
        return self

    # ------------------------------------------------------------------------- anchor
    def anchor_terms(self, X: torch.Tensor):
        return self.anchor(X, self.candidates)

    def additive_coefficients(self, X: torch.Tensor) -> dict:
        """{edge: float} for the support oracle. Fixed for this X, whatever the mask."""
        with torch.no_grad():
            _, a = self.anchor_terms(X)
        return {e: float(a[k]) for k, e in enumerate(self.candidates)}

    def anchor_bias(self, X: torch.Tensor) -> float:
        with torch.no_grad():
            b, _ = self.anchor_terms(X)
        return float(b)

    # ------------------------------------------------------------------------- P and A
    def indicator(self, structure: Iterable[tuple]) -> torch.Tensor:
        z = torch.zeros(len(self.candidates), dtype=self.dtype)
        for e in structure:
            e = tuple(e)
            if e not in self.index:
                raise KeyError(f"edge {e} is not in the candidate universe")
            z[self.index[e]] = 1.0
        return z

    def adjacency(self, z: torch.Tensor, base: torch.Tensor) -> torch.Tensor:
        """B + P z, symmetric, differentiable in z. No normalisation is applied."""
        if len(self.candidates) == 0:
            return base
        delta = torch.zeros_like(base)
        delta = delta.index_put((self._i, self._j), z, accumulate=True)
        delta = delta.index_put((self._j, self._i), z, accumulate=True)
        return base + delta

    # ------------------------------------------------------------------------ forward
    def linear_value(self, X: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        b, a = self.anchor_terms(X)
        return b + (a * z).sum()

    def forward(self, X: torch.Tensor, z: torch.Tensor,
                base: torch.Tensor | None = None) -> torch.Tensor:
        lin = self.linear_value(X, z)
        if not self.residual_active:
            return lin                       # EXACT bypass: g is never evaluated
        if base is None:
            base = backbone_adjacency(self.n_nodes, dtype=self.dtype)
        return lin + self.gamma * self.g(X, self.adjacency(z, base))

    # --------------------------------------------------------------------- diagnostics
    def interaction_gap(self, X, z_a, z_b, base=None) -> float:
        """Second-order interaction  f(a+b) - f(a) - f(b) + f(0).

        Exactly 0 for any additive f, so a non-zero value is behavioural evidence that
        the residual is actually exercised — unlike a weight norm, which says nothing
        about whether the nonlinearity ever changes an output.

        REQUIRES DISJOINT SUPPORTS. If z_a and z_b share an edge then the union is not
        their sum, and the identity fails even for a purely additive f; the quantity
        would then measure the overlap, not the interaction.
        """
        if float(torch.dot(z_a, z_b)) != 0.0:
            raise ValueError("interaction_gap requires disjoint supports; "
                             "z_a and z_b share at least one selected edge")
        z0 = torch.zeros_like(z_a)
        z_ab = z_a + z_b
        with torch.no_grad():
            return float(self(X, z_ab, base) - self(X, z_a, base)
                         - self(X, z_b, base) + self(X, z0, base))
