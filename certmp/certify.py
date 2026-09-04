"""Certificates a designer would act on."""
from .reach import exact_interval

def certify_threshold(model, n, mandatory, optional, X, tau):
    """Sound certificate: no structure in the uncertainty set pushes the score above tau."""
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
