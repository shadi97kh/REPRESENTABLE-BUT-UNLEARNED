"""F4: exact deterministic top-k rank stability for a candidate library."""
import random
from certmp.models import MonoMPNN
from certmp.ensemble import bpp_lattice, onehot_features
from certmp.reach import exact_interval
from certmp.certify import certify_topk
from certmp.provenance import new_run, record

def main():
    run = new_run("f4_topk_stability", dict(library=40, length=50, k=5))
    random.seed(2)
    m = None; cands = []
    for t in range(40):
        seq = "".join(random.choice("ACGU") for _ in range(50))
        lat = bpp_lattice(seq)
        if lat["k"] == 0: continue
        X = onehot_features(seq)
        if m is None: m = MonoMPNN(X.shape[1], 8, 2, agg="sum", nonneg=True, seed=0)
        lo, hi, _ = exact_interval(m, lat["n"], lat["mandatory"], lat["optional"], X)
        cands.append({"name": f"cand{t:02d}", "lo": lo, "hi": hi, "k": lat["k"]})
    for K in (1, 3, 5, 10):
        r = certify_topk(cands, K)
        print(f"top-{K:<3d} certified={str(r['certified']):5s} margin={r['margin']:12.4f}  {r['topk'][:3]}")
    record(run, "f4_topk_stability",
           {"candidates": cands, "results": {K: certify_topk(cands, K) for K in (1,3,5,10)}})

if __name__ == "__main__": main()
