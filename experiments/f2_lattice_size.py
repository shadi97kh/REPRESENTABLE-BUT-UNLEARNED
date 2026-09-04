"""F2: registered prediction R7. Does the uncertainty lattice have practical bite?
If k is small on real ensembles the theorem is true but vacuous."""
import random, statistics
from certmp.ensemble import bpp_lattice
from certmp.provenance import new_run, record

LENGTHS = (40, 80, 150)
N_SEQ = 30

def main():
    run = new_run("f2_lattice_size", dict(lengths=LENGTHS, n_seq=N_SEQ, lo=0.05, hi=0.90))
    random.seed(0); out = {}
    for L in LENGTHS:
        ks = []
        for _ in range(N_SEQ):
            s = "".join(random.choice("ACGU") for _ in range(L))
            ks.append(bpp_lattice(s)["k"])
        ks.sort()
        out[L] = {"median_k": statistics.median(ks), "p90_k": ks[int(.9*len(ks))],
                  "max_k": max(ks), "median_lattice_log10": statistics.median(ks)*0.30103}
        print(f"len={L:4d}  median k={statistics.median(ks):5.1f}  p90={ks[int(.9*len(ks))]:4d}  "
              f"max={max(ks):4d}  median 2^k = 10^{statistics.median(ks)*0.30103:.1f}")
    vacuous = out[150]["median_k"] < 8
    print(f"\nR7 {'FAILS -- result is vacuous' if vacuous else 'holds -- brute force intractable'}")
    record(run, "f2_lattice_size", {"by_length": out, "vacuous": vacuous})

if __name__ == "__main__": main()
