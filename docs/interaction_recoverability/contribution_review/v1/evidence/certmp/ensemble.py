"""Real RNA structural uncertainty from the Boltzmann ensemble via ViennaRNA.

Current siRNA pipelines take the MFE structure and treat it as certain. The ensemble
says otherwise. Base pairs are split by their equilibrium probability:
    p > hi     -> mandatory  (present in essentially every structure)
    lo <= p <= hi -> optional (the uncertainty lattice)
    p < lo     -> discarded
BANDED LATTICE, NOT SOUND. This function keeps its original probability band for the
experiments that were run against it (f2, f3, f3b, f6). It is NOT sound over the ensemble:
it forces mandatory pairs present, which invalidates any lower bound, and it discards pairs
below `lo`, which leaves structures outside the lattice and therefore outside the upper
bound. At lo=0.05, f7 measured 73.8% coverage and observed 9 actual upper-bound violations
in 20000 sampled structures.

USE sound_lattice() BELOW for any claim about the Boltzmann ensemble.
"""
def bpp_lattice(seq, lo=0.05, hi=0.90, floor=1e-3):
    import RNA
    fc = RNA.fold_compound(seq)
    ss, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    mandatory, optional, probs = [], [], {}
    for i in range(1, n+1):
        for j in range(i+1, n+1):
            p = bpp[i][j]
            if p < floor: continue
            e = (i-1, j-1)                      # 0-indexed nodes
            probs[e] = p
            if p > hi: mandatory.append(e)
            elif p >= lo: optional.append(e)
    return {"n": n, "mfe_structure": ss, "mfe_energy": mfe,
            "mandatory": mandatory, "optional": optional, "probs": probs,
            "k": len(optional), "lattice_size_log10": len(optional)*0.30103}

CANONICAL = {("A", "U"), ("U", "A"), ("G", "C"), ("C", "G"), ("G", "U"), ("U", "G")}
MIN_LOOP = 3          # ViennaRNA TURN: a hairpin needs at least 3 unpaired bases


def sound_lattice(seq, floor=0.0, canonical_only=True):
    """SOUND lattice: no mandatory edges, and no pair discarded.

    This is the fix for the two soundness defects f3b measured.
      - mandatory = [] means nothing is forced present, so the empty-optional endpoint is
        a subset of every real structure's edge set.
      - floor = 0.0 keeps EVERY pair, so the lattice is the full power set of base pairs
        and every secondary structure over this sequence is inside it by construction.
    With both, monotonicity gives model(structure) <= model(all edges) for every structure
    in the ensemble, with no coverage caveat. floor > 0 trades that guarantee for a smaller
    k and a tighter bound; experiments/f7_soundness.py measures the trade.

    canonical_only=True keeps soundness while shrinking the lattice for free. ViennaRNA's
    energy model can only form the six canonical pairs, and only with at least MIN_LOOP
    unpaired bases between partners, so no structure in the ensemble can contain any pair
    this excludes. Excluding them therefore removes edges no real structure uses rather
    than edges the bound needs. Set it False to get the literal full power set.
    """
    import RNA
    fc = RNA.fold_compound(seq)
    ss, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    bpp = fc.bpp()
    n = len(seq)
    u = seq.replace("T", "U").upper()
    optional, probs = [], {}
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            if canonical_only:
                if j - i <= MIN_LOOP: continue
                if (u[i - 1], u[j - 1]) not in CANONICAL: continue
            p = bpp[i][j]
            if p < floor: continue          # floor=0.0 keeps every pair: p >= 0 always
            e = (i - 1, j - 1)
            probs[e] = p
            optional.append(e)
    return {"n": n, "mfe_structure": ss, "mfe_energy": mfe,
            "mandatory": [], "optional": optional, "probs": probs,
            "k": len(optional), "floor": floor, "canonical_only": canonical_only,
            "lattice_size_log10": len(optional) * 0.30103,
            "sound_over_ensemble": floor <= 0.0}


def generic_features(seq, decay=3.0):
    """Non-negative, LENGTH-INDEPENDENT node features (H4).

    f6's joint (position, nucleotide) one-hot cannot transfer between lengths, but the
    target-context experiment needs one model applied to 21, 50, 100 and 150 nt windows.
    This keeps the dominant positional effect, end asymmetry, in a length-free form:
    nucleotide identity, and nucleotide weighted by proximity to each end. All entries
    lie in [0, 1], so H4 holds.
    """
    import numpy as np
    idx = {c: t for t, c in enumerate("ACGU")}
    L = len(seq); X = np.zeros((L, 12))
    for t, c in enumerate(seq.replace("T", "U")):
        if c not in idx: continue
        j = idx[c]
        w5 = float(np.exp(-t / decay))
        w3 = float(np.exp(-(L - 1 - t) / decay))
        X[t, j] = 1.0
        X[t, 4 + j] = w5
        X[t, 8 + j] = w3
    return X


def onehot_features(seq, extra=None):
    """Non-negative node features (hypothesis H4). One-hot nucleotide + optional
    non-negative scalars such as unpaired probability."""
    import numpy as np
    idx = {c: t for t, c in enumerate("ACGU")}
    X = np.zeros((len(seq), 4))
    for t, c in enumerate(seq.replace("T", "U")):
        if c in idx: X[t, idx[c]] = 1.0
    if extra is not None:
        E = np.asarray(extra, dtype=float).reshape(len(seq), -1)
        assert E.min() >= 0, "extra features must be non-negative (H4)"
        X = np.concatenate([X, E], axis=1)
    return X
