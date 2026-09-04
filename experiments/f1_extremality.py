"""F1: the theorem. Registered predictions R1-R6.
Brute-force oracle over 2^k states vs the 2-pass endpoint bound."""
import numpy as np
from certmp.models import MonoMPNN
from certmp.reach import endpoint_gap
from certmp.provenance import new_run, record

def trial(agg, nonneg, neg_X, seed, n=7, k=8, d_in=3):
    r = np.random.RandomState(seed)
    mandatory = [(0,1),(1,2),(2,3)]
    cand = [(i,j) for i in range(n) for j in range(i+1,n) if (i,j) not in mandatory]
    optional = [cand[t] for t in r.choice(len(cand), size=k, replace=False)]
    X = r.randn(n, d_in) if neg_X else r.rand(n, d_in)*2.0
    m = MonoMPNN(d_in, 8, 2, agg=agg, nonneg=nonneg, seed=seed)
    return endpoint_gap(m, n, mandatory, optional, X)

CASES = [("R1 sum   nonnegW nonnegX", "sum",     True,  False, True),
         ("R2 max   nonnegW nonnegX", "max",     True,  False, True),
         ("R3 mean  nonnegW nonnegX", "mean",    True,  False, False),
         ("R4 degnorm nonnegW nonnegX","degnorm",True,  False, False),
         ("R5 sum   SIGNEDW nonnegX", "sum",     False, False, False),
         ("R6 sum   nonnegW NEG X",   "sum",     True,  True,  False)]
TRIALS = 25

def main():
    run = new_run("f1_extremality", dict(trials=TRIALS, n=7, k=8))
    out = []
    print(f"{'case':30s} {'exact':>10s} {'max gap_hi':>13s} {'max gap_lo':>13s}  verdict")
    for label, agg, nonneg, negX, expect in CASES:
        rs = [trial(agg, nonneg, negX, s) for s in range(TRIALS)]
        ne = sum(r["exact"] for r in rs)
        ghi = max(r["gap_hi"] for r in rs); glo = max(r["gap_lo"] for r in rs)
        ok = (ne == TRIALS) == expect
        print(f"{label:30s} {ne:4d}/{TRIALS:<5d} {ghi:13.3e} {glo:13.3e}  {'PASS' if ok else '*** FAIL ***'}")
        out.append(dict(case=label, agg=agg, nonneg=nonneg, neg_X=negX, n_exact=ne,
                        trials=TRIALS, max_gap_hi=ghi, max_gap_lo=glo,
                        expected_exact=expect, passed=ok, states_per_trial=rs[0]["n_states"]))
    record(run, "f1_extremality", {"cases": out,
        "theorem_holds": all(o["passed"] for o in out)})

if __name__ == "__main__": main()
