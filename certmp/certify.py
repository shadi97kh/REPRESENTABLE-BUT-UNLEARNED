"""Certificates a designer would act on."""
from .reach import exact_interval

def certify_threshold(model, n, mandatory, optional, X, tau):
    """Certificate: no edge subset in the LATTICE pushes the score above tau.

    Scope, measured in experiments/f3b_relaxation_slack.py. 'worst_case' is a sound upper
    bound over the lattice, and over every real structure contained in the lattice, but it
    is loose: a median 1.82x the true sampled Boltzmann maximum. Structures using pairs
    dropped below the `lo` band are not covered unless the lattice was built with lo=0.

    'best_case' is NOT a sound lower bound on the ensemble. It forces all mandatory pairs
    present, and real structures omit them. Reported for lattice diagnostics only.
    """
    lo, hi, passes = exact_interval(model, n, mandatory, optional, X)
    return {"certified": hi <= tau, "worst_case": hi, "best_case": lo, "tau": tau,
            "forward_passes": passes, "lattice_size_log10": len(optional)*0.30103,
            "margin": tau - hi}

def certify_topk(candidates, k):
    """Exact, DETERMINISTIC top-k stability. candidates: list of dicts with
    'name', 'lo', 'hi' from exact_interval. The ranking is certified when the k-th
    best worst-case strictly exceeds the (k+1)-th best best-case: intervals do not
    straddle the cut, so no perturbation in the set can reorder across it.
    Contrast: Jia et al. ICLR 2020 give only a probabilistic smoothing certificate."""
    s = sorted(candidates, key=lambda c: -c["lo"])
    if len(s) <= k: return {"certified": True, "topk": [c["name"] for c in s], "margin": float("inf")}
    margin = s[k-1]["lo"] - s[k]["hi"]
    return {"certified": margin > 0, "topk": [c["name"] for c in s[:k]],
            "margin": margin, "cut_lo": s[k-1]["lo"], "cut_hi": s[k]["hi"]}
