"""Adapters that re-index an oracle's ground set without changing its family.

`CertifiedModel` indexes candidates as node pairs `(i, j)`. Some families name their
ground set differently -- a layered DAG names edges `(layer, u, v)`. Re-labelling is a
bijection on the ground set, so the family, its counts and its supports are unchanged;
only the names move. This is the adapter the pilot dry run identified as a blocker.
"""
from __future__ import annotations

from .base import Feasibility, OracleResult, SupportOracle, normalise_mask


class EdgeRelabeledOracle(SupportOracle):
    """Wrap `base` so its ground set is presented under new names.

    `mapping` must be a bijection from the base ground set onto the new names; the
    constructor checks that, because a collision would silently merge two distinct edges
    and change the family.
    """

    def __init__(self, base: SupportOracle, mapping: dict):
        self.base = base
        old = tuple(base.ground_set())
        if set(mapping) != set(old):
            raise ValueError("mapping must cover exactly the base ground set")
        new = [tuple(mapping[e]) for e in old]
        if len(set(new)) != len(new):
            raise ValueError("mapping is not injective: two edges share a new name")
        self.fwd = {e: tuple(mapping[e]) for e in old}
        self.bwd = {tuple(mapping[e]): e for e in old}
        self._ground = tuple(new)
        self.family_id = f"{base.family_id}|relabeled"

    def _to_base(self, edges):
        return tuple(self.bwd[tuple(e)] for e in edges)

    def _from_base(self, edges):
        return tuple(sorted(self.fwd[tuple(e)] for e in edges))

    def ground_set(self):
        return self._ground

    def assumptions(self):
        return self.base.assumptions() + ("ground set re-labelled by a bijection",)

    def describe(self):
        d = dict(self.base.describe())
        d["relabeled_to"] = [list(e) for e in self._ground]
        return d

    def feasibility(self, forced=(), forbidden=()):
        f, b = normalise_mask(forced, forbidden)
        return self.base.feasibility(self._to_base(f), self._to_base(b))

    def count(self, forced=(), forbidden=()):
        f, b = normalise_mask(forced, forbidden)
        return self.base.count(self._to_base(f), self._to_base(b))

    def log_partition(self, coefficients, forced=(), forbidden=()):
        f, b = normalise_mask(forced, forbidden)
        c = {self.bwd[tuple(k)]: v for k, v in (coefficients or {}).items()}
        return self.base.log_partition(c, self._to_base(f), self._to_base(b))

    def support(self, coefficients, forced=(), forbidden=()):
        f, b = normalise_mask(forced, forbidden)
        c = {self.bwd[tuple(k)]: v for k, v in (coefficients or {}).items()}
        r = self.base.support(c, self._to_base(f), self._to_base(b))
        if r.witness is None:
            return r
        return OracleResult(r.value, self._from_base(r.witness), r.numerical_status,
                            r.termination_status, r.family_hash, r.bound_direction,
                            self.assumptions(), r.stats)

    def stats(self):
        return self.base.stats()
