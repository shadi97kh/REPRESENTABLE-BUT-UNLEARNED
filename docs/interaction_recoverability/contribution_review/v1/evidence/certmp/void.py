class VoidRun(Exception): pass

def check_finite(*arrays):
    import numpy as np
    for a in arrays:
        if not np.all(np.isfinite(a)): raise VoidRun("non-finite value")

def check_live(values, tol=1e-12):
    """Output must vary across the uncertainty set, else 'exactness' is vacuous."""
    lo, hi = min(values), max(values)
    if hi - lo <= tol: raise VoidRun(f"dead network: output constant across set ({lo:.3e})")

def check_nonempty_lattice(optional):
    if len(optional) == 0: raise VoidRun("k=0: empty uncertainty lattice, nothing to certify")

def relgap(a, b, scale):
    return (a - b) / max(abs(scale), 1e-12)
