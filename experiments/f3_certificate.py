"""F3: end-to-end certificate on a real folded sequence.

RETIRED, kept for provenance. The model is untrained, so every magnitude here is an
initialisation artifact; f6 replaces them. The lattice is the BANDED one, which f7 shows
is not sound over the ensemble; f7 replaces it. Do not quote this experiment.

Reports the MFE-only prediction alongside the certified worst case, i.e. exactly what
current pipelines miss by treating the MFE structure as certain."""
import random, numpy as np
from certmp.models import MonoMPNN
from certmp.ensemble import bpp_lattice, onehot_features
from certmp.certify import certify_threshold
from certmp.reach import build_A
from certmp.provenance import new_run, record

def main():
    run = new_run("f3_certificate", dict(lo=0.05, hi=0.90, length=60, n_seq=20))
    random.seed(1); rows = []
    for t in range(20):
        seq = "".join(random.choice("ACGU") for _ in range(60))
        lat = bpp_lattice(seq)
        if lat["k"] == 0: continue
        X = onehot_features(seq)
        m = MonoMPNN(X.shape[1], 8, 2, agg="sum", nonneg=True, seed=t)
        # MFE-only prediction: the single minimum-free-energy structure
        mfe_edges = []
        stack = []
        for i, c in enumerate(lat["mfe_structure"]):
            if c == "(": stack.append(i)
            elif c == ")": mfe_edges.append((stack.pop(), i))
        y_mfe = m.forward(X, build_A(lat["n"], mfe_edges, [], []))
        cert = certify_threshold(m, lat["n"], lat["mandatory"], lat["optional"], X, tau=1e9)
        under = (cert["worst_case"] - y_mfe) / max(abs(y_mfe), 1e-12)
        rows.append(dict(k=lat["k"], mfe_energy=lat["mfe_energy"], y_mfe=y_mfe,
                         worst_case=cert["worst_case"],
                         mfe_underestimate_rel=under, passes=cert["forward_passes"],
                         lattice_log10=cert["lattice_size_log10"]))
        print(f"seq{t:02d} k={lat['k']:3d} 2^k=10^{cert['lattice_size_log10']:5.1f}  "
              f"MFE={y_mfe:10.3f}  certified worst={cert['worst_case']:10.3f}  "
              f"MFE underestimates by {under*100:6.1f}%  passes={cert['forward_passes']}")
    med = float(np.median([r["mfe_underestimate_rel"] for r in rows])) if rows else float("nan")
    print(f"\nmedian relative underestimate of worst-case risk by MFE-only: {med*100:.1f}%")
    record(run, "f3_certificate", {"rows": rows, "median_mfe_underestimate_rel": med})

if __name__ == "__main__": main()
