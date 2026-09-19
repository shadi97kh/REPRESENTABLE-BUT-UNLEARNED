"""A toy target with genuinely non-additive interactions over a FEASIBLE family.

The point is to have a target that no additive model can fit, where the interaction is
between edges that can actually co-occur in a valid structure. A target built from
crossing pairs would be vacuous: those structures are not in the family, so an additive
model would never be penalised for missing them.

STACKING. Base pairs (i, j) and (i+1, j-1) stack. Stacking is the dominant stabilising
interaction in real RNA thermodynamics, it is a second-order interaction between two
edges, and the two edges are nested, so every stacked pair is feasible.

    y(z) = sum_e w_e z_e  +  lambda * #{ (i,j) in z : (i+1,j-1) in z }

The first term is additive and any additive model can fit it. The second cannot be
written as sum_e c_e z_e, which `additive_least_squares_floor` demonstrates by fitting
the *unconstrained* best additive model and reporting its residual. That floor is the
baseline to beat, and it is the strongest additive baseline available — not a weakened
one.
"""
from __future__ import annotations

import random

import numpy as np

from .oracles import RNANonCrossingOracle
from .oracles.bruteforce import enumerate_rna


def toy_features(seq: str):
    """One-hot nucleotide concatenated with one-hot position.

    Nucleotide identity alone is NOT enough: the anchor computes a_e from the endpoint
    attributes, so with a 4-dim one-hot every candidate sharing an endpoint nucleotide
    pair is forced to the same coefficient. On a sequence like GGGGAAAACCCC every
    candidate is a G-C pair, which collapses the whole anchor to a single number and
    makes any comparison against a per-edge additive baseline meaningless. Position
    restores the ability to express an arbitrary a_e, which is the encoding cost of
    keeping the anchor exactly additive.
    """
    import torch
    n = len(seq)
    X = torch.zeros(n, 4 + n, dtype=torch.float64)
    for t, c in enumerate(seq):
        X[t, "ACGU".index(c)] = 1.0
        X[t, 4 + t] = 1.0
    return X


def count_stacks(structure) -> int:
    S = {tuple(e) for e in structure}
    return sum(1 for (i, j) in S if (i + 1, j - 1) in S)


def stacking_toy(seq: str = "GGGGAAAACCCC", lam: float = 3.0, seed: int = 0,
                 min_loop: int = 3, canonical_only: bool = True):
    """(oracle, candidates, members, targets, meta). Targets are exact, not sampled."""
    o = RNANonCrossingOracle(seq, min_loop, canonical_only)
    cand = o.ground_set()
    members = [tuple(sorted(s)) for s in enumerate_rna(seq, min_loop, canonical_only)]
    rng = random.Random(seed)
    w = {e: rng.uniform(-1.0, 1.0) for e in cand}
    targets = {s: sum(w[e] for e in s) + lam * count_stacks(s) for s in members}
    meta = {"seq": seq, "lam": lam, "seed": seed, "n_members": len(members),
            "n_candidates": len(cand),
            "stack_counts": sorted({count_stacks(s) for s in members})}
    return o, cand, members, targets, meta


def design_matrix(structures, candidates) -> np.ndarray:
    """[1, z] — the feature map of the most general additive model."""
    return np.array([[1.0] + [1.0 if e in set(s) else 0.0 for e in candidates]
                     for s in structures])


def additive_least_squares_floor(train, val, candidates, targets):
    """Best *unconstrained* additive model, fitted on train, scored on both splits.

    This is the strongest additive baseline: a free parameter per candidate edge plus an
    intercept, fitted in closed form. Any additive model, however parameterised, has
    train MSE at least this large.
    """
    Dtr, Dva = design_matrix(train, candidates), design_matrix(val, candidates)
    ytr = np.array([targets[s] for s in train])
    yva = np.array([targets[s] for s in val])
    coef, *_ = np.linalg.lstsq(Dtr, ytr, rcond=None)
    return {"train_mse": float(((Dtr @ coef - ytr) ** 2).mean()),
            "val_mse": float(((Dva @ coef - yva) ** 2).mean()),
            "coef": coef}


def split(members, frac_train: float = 0.7, seed: int = 11):
    idx = list(range(len(members)))
    random.Random(seed).shuffle(idx)
    k = int(frac_train * len(members))
    return [members[i] for i in idx[:k]], [members[i] for i in idx[k:]]
