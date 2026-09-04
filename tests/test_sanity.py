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

def test_torch_numpy_parity():
    """The certificate is computed on the numpy model. If the torch model it was
    trained as disagrees, the certificate describes a different network."""
    try:
        import torch
    except ImportError:
        print("  (torch absent, parity test skipped)"); return
    from certmp.models import TorchMonoMPNN
    rs = np.random.RandomState(3)
    for agg in ("sum", "max"):
        for layers in (1, 2, 3):
            tm = TorchMonoMPNN(4, 7, layers, agg=agg, seed=layers)
            nm = tm.export_numpy()
            n = 9
            X = rs.rand(n, 4)
            A = np.eye(n)
            for _ in range(12):
                i, j = rs.randint(0, n, 2); A[i, j] = A[j, i] = 1.0
            with torch.no_grad():
                y_t = float(tm(torch.tensor(X, dtype=torch.float32).unsqueeze(0),
                               torch.tensor(A, dtype=torch.float32).unsqueeze(0))[0])
            y_n = nm.forward(X, A)
            rel = abs(y_t - y_n) / max(abs(y_n), 1e-12)
            assert rel < 1e-5, f"{agg}/{layers}L parity broken: torch {y_t} vs numpy {y_n} (rel {rel:.2e})"
            # a NON-TRIVIAL affine readout must also survive the export
            tm.set_output_affine(3.7e-4, -2.5)
            nm2 = tm.export_numpy()
            with torch.no_grad():
                y_t2 = float(tm(torch.tensor(X, dtype=torch.float32).unsqueeze(0),
                                torch.tensor(A, dtype=torch.float32).unsqueeze(0))[0])
            y_n2 = nm2.forward(X, A)
            rel2 = abs(y_t2 - y_n2) / max(abs(y_n2), 1e-12)
            assert rel2 < 1e-5, f"{agg}/{layers}L affine parity broken ({rel2:.2e})"

def test_affine_readout_preserves_extremality():
    """scale >= 0 is a monotone map and shift is order-preserving, so the endpoint
    theorem must be unaffected. If it is not, the affine readout broke the certificate."""
    from certmp.reach import endpoint_gap
    rs = np.random.RandomState(11)
    n, k = 6, 7
    mand = [(0,1),(1,2)]
    opt = [(0,3),(1,4),(2,5),(3,4),(4,5),(0,5),(2,3)][:k]
    X = rs.rand(n, 3) * 2
    for agg in ("sum", "max"):
        for scale, shift in ((1e-4, -3.0), (5.0, 2.0), (0.0, 1.5)):
            m = MonoMPNN(3, 8, 2, agg=agg, nonneg=True, seed=1)
            m.scale, m.shift = scale, shift
            if scale == 0.0:
                continue                      # constant output would VOID the trial
            g = endpoint_gap(m, n, mand, opt, X)
            assert g["exact"], f"{agg} scale={scale} shift={shift}: extremality lost {g}"

def test_negative_scale_rejected():
    m = MonoMPNN(3, 8, 2, agg="sum", nonneg=True, seed=0)
    m.scale = -1.0
    assert not m.structurally_certifiable(), "negative output scale must leave the safe class"

if __name__ == "__main__":
    test_exact_matches_brute_force(); test_rejects_uncertified_class()
    test_rejects_negative_features(); test_void_on_empty_lattice()
    test_torch_numpy_parity()
    test_affine_readout_preserves_extremality(); test_negative_scale_rejected()
    print("sanity: all passed")
