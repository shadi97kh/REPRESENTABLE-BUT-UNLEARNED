"""Necessary deterministic checks of new conventions and public schedules.

No random generation, sampled-data estimator, fit, or previous check suite runs.
The four hand-written ordinates below test interpolation arithmetic only.
"""
import importlib.util
import json
import os
from pathlib import Path
import sys
import time
import mpmath as mp

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE))
import estimator as e


def main():
    start = time.process_time()
    assert len(os.sched_getaffinity(0)) == 1
    assert os.environ.get("CUDA_VISIBLE_DEVICES") == ""
    for key in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
        assert os.environ.get(key) == "1"
    mp.mp.dps = 70
    previous = json.loads((ROOT / "docs/interaction_recoverability/local_five_law/v1/checks.json").read_text())
    assert previous["status"] == "passed"
    q = e.linear_quantile_from_sorted([1,4,9,16])
    points = ["0", ".125", ".25", ".375", ".5", ".875", "1"]
    expected = ["1", "1", "1", "2.5", "4", "12.5", "16"]
    actual = [q(mp.mpf(p)) for p in points]
    assert actual == [mp.mpf(v) for v in expected]
    try:
        e.estimate_from_sources(object())
    except PermissionError:
        guard = "blocked before reading source input"
    else:
        raise AssertionError("Sample execution guard failed")
    alpha = [mp.mpf(2)/3, mp.mpf(3)/4]
    beta = [mp.mpf(3)/5, mp.mpf(4)/5]
    aa, bb = [0,*alpha,1,0], [0,*beta,0,1]
    curves = [(lambda z, a=a, b=b: mp.exp(z+a)+mp.exp(e.h(z)+b))
              for a,b in zip(aa,bb)]
    op = e.AnchorOperators(curves, depth=8, tolerance=mp.mpf("1e-40"))
    g10, g20 = op.g1(0), op.g2(0)
    assert abs(g10-1) < mp.mpf("1e-25")
    assert abs(g20-1) < mp.mpf("1e-25")
    # This checks the population remainder's sign and its normalization,
    # not a new convergence scan or the previous tangent-check suite.
    x, y = mp.mpf("4.5"), mp.mpf(4)
    xk, yk = x, y
    for _ in range(op.depth):
        xk, yk = op.contract(xk), op.contract(yk)
    finite = op.difference(x,y)
    exact_finite = mp.exp(x)-mp.exp(y)-(mp.exp(xk)-mp.exp(yk))
    remainder_error = abs(finite-exact_finite)
    assert remainder_error < mp.mpf("1e-25")
    domain = []
    for n in [10**12, 10**32, 10**64]:
        s = e.schedule(5*n)
        allowed, rank, required = e.domain_condition(s)
        assert mp.power(mp.mpf(575)/613,s.depth) <= mp.mpf(s.N)**-2
        domain.append(dict(n=str(n), N=str(s.N), depth=s.depth,
                           lattice=s.lattice, radius=mp.nstr(s.radius,16),
                           lower_effective_rank=mp.nstr(rank,16),
                           required_rank=mp.nstr(required,16),
                           domain_condition=bool(allowed),
                           precision_bits=s.precision_bits,
                           quadrature_mesh=mp.nstr(s.mesh,8)))
    out = dict(
        status="passed", backend="mpmath, 70 decimal digits; diagnostic, not certified interval arithmetic",
        previous_checks="read existing passed record only; not rerun",
        interpolation=dict(probabilities=points, values=[str(x) for x in actual]),
        sample_entry=guard,
        population_normalization=dict(g1_zero=mp.nstr(g10,30), g2_zero=mp.nstr(g20,30)),
        finite_population_remainder_error=mp.nstr(remainder_error,20),
        domain_checks=domain,
        scope="Only deterministic interpolation, normalization, finite-population identity and public-count feasibility checks. No sample estimator, data generation, variance or performance experiment.",
        cpu_s=time.process_time()-start, numerical_workers=1, numerical_threads=1,
        affinity=sorted(os.sched_getaffinity(0)), gpu_s=0)
    with (HERE/"checks.json").open("x") as f:
        json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))


if __name__ == "__main__":
    main()
