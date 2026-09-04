import numpy as np
from certmp.models import MonoMPNN
from certmp.reach import exact_interval, brute_force_interval, build_A
from certmp.void import VoidRun

def test_exact_matches_brute_force():
    n, k = 6, 7
    mand = [(0,1),(1,2)]
    opt = [(0,3),(1,4),(2,5),(3,4),(4,5),(0,5),(2,3)][:k]
    X = np.random.RandomState(0).rand(n, 3)*2
    for agg in ("sum","max"):
        m = MonoMPNN(3, 8, 2, agg=agg, nonneg=True, seed=0)
        lo, hi, p = exact_interval(m, n, mand, opt, X)
        blo, bhi, ns = brute_force_interval(m, n, mand, opt, X)
        assert abs(lo-blo) < 1e-9 and abs(hi-bhi) < 1e-9, f"{agg}: 2-pass != brute force"
        assert p == 2 and ns == 2**k

def test_rejects_uncertified_class():
    X = np.random.RandomState(0).rand(5,3)
    for agg, nonneg in (("mean",True),("degnorm",True),("sum",False)):
        m = MonoMPNN(3, 8, 2, agg=agg, nonneg=nonneg, seed=0)
        try:
            exact_interval(m, 5, [(0,1)], [(1,2)], X); assert False, f"{agg}/{nonneg} not rejected"
        except ValueError: pass

def test_rejects_negative_features():
    m = MonoMPNN(3, 8, 2, agg="sum", nonneg=True, seed=0)
    X = -np.ones((5,3))
    try:
        exact_interval(m, 5, [(0,1)], [(1,2)], X); assert False, "negative X not rejected"
    except ValueError: pass

def test_void_on_empty_lattice():
    m = MonoMPNN(3, 8, 2, agg="sum", nonneg=True, seed=0)
    X = np.random.RandomState(0).rand(5,3)
    try:
        brute_force_interval(m, 5, [(0,1)], [], X); assert False, "empty lattice not flagged"
    except VoidRun: pass

if __name__ == "__main__":
    test_exact_matches_brute_force(); test_rejects_uncertified_class()
    test_rejects_negative_features(); test_void_on_empty_lattice()
    print("sanity: all passed")
