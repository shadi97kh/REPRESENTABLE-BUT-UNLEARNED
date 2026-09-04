"""Real RNA structural uncertainty from the Boltzmann ensemble via ViennaRNA.

Current siRNA pipelines take the MFE structure and treat it as certain. The ensemble
says otherwise. Base pairs are split by their equilibrium probability:
    p > hi     -> mandatory  (present in essentially every structure)
    lo <= p <= hi -> optional (the uncertainty lattice)
    p < lo     -> discarded
The lattice is a SOUND OUTER RELAXATION of the ensemble: not every edge subset is a valid
secondary structure (nesting, one-partner constraints), so certifying over the lattice
certifies over the ensemble, conservatively. Quantify that conservatism, never hide it.
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
