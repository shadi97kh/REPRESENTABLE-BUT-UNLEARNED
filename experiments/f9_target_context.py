"""F9 (step 3): run the SOUND certificate on real mRNA target-site windows.

f6 exposed a scope problem. An isolated 21 nt siRNA duplex has a median of 6 optional
pairs, so its uncertainty set holds about 64 structures and brute force is affordable.
The exponential-to-constant collapse only earns its keep where the uncertainty set is
large, which f2 located in longer sequence: k reaches 86 at 150 nt.

This runs the certificate where it belongs, on the mRNA target site the siRNA binds,
extracted from GENCODE transcripts for the Huesken genes, at 50, 100 and 150 nt.

Lattice: mandatory = the backbone, which every structure has by construction and which the
model saw in training; optional = every canonical base pair with the minimum hairpin loop,
none discarded. No base pair is forced, so the upper bound is sound over the ensemble, and
coverage against pbacktrack samples is checked rather than assumed.

Model: the length-generalizing monotone network (data/models/monotone_generic.json,
held-out Spearman 0.385). f8 shows this encoding costs about 0.21 Spearman against an
unconstrained model, so treat magnitudes here as indicative of the certificate's behaviour,
not as a state-of-the-art efficacy predictor.
"""
import gzip, os, statistics
import numpy as np
import RNA
from certmp.ensemble import sound_lattice, generic_features
from certmp.reach import exact_interval, build_A
from certmp.train import load_numpy_model
from certmp.provenance import new_run, record
from experiments._huesken import pairs_of, load_verified, gene_map

GENCODE = os.environ.get("GENCODE_FA", "data/gencode/gencode.v47.pc_transcripts.fa.gz")
MODEL_PATH = "data/models/monotone_generic.json"
LENGTHS = (50, 100, 150)
N_SITES = 24
N_SAMPLES = 300
RNG_SEED = 20260904

# Huesken-era symbols that GENCODE now files under a current name. Rat genes (lowercase
# in the source annotation) are not in a human build and are expected to be unmatched.
ALIASES = {"HIP2": "UBE2K", "HSPC150": "UBE2T", "RAB6IP1": "DENND5A", "TC10": "RHOQ",
           "UBE2N": "UBE2N", "P2RX3": "P2RX3"}


def rc(s):
    return s.translate(str.maketrans("ACGU", "UGCA"))[::-1]


def load_transcripts(path, wanted):
    """GENCODE fasta headers are pipe-delimited with the gene symbol in field 6."""
    want = {w.upper() for w in wanted}
    out, name, buf = {}, None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if name and name in want:
                    out.setdefault(name, []).append("".join(buf).upper().replace("T", "U"))
                f = line[1:].split("|")
                name = f[5].upper() if len(f) > 5 else None
                buf = []
            elif name in want:
                buf.append(line.strip())
        if name and name in want:
            out.setdefault(name, []).append("".join(buf).upper().replace("T", "U"))
    return out


def sample_max(seq, m, X, backbone, n, n_samples):
    md = RNA.md(); md.uniq_ML = 1
    fc = RNA.fold_compound(seq, md)
    _, e = fc.mfe(); fc.exp_params_rescale(e); fc.pf()
    ys, sets = [], []
    for s in fc.pbacktrack(n_samples):
        P = pairs_of(s); sets.append(P)
        ys.append(m.forward(X, build_A(n, backbone + sorted(P), [], [])))
    return float(max(ys)), sets


def main():
    for p in (GENCODE, MODEL_PATH):
        if not os.path.exists(p):
            print(f"MISSING: {p}. Run `make gencode` / `make f6`."); raise SystemExit(2)
    run = new_run("f9_target_context", dict(gencode=GENCODE, lengths=LENGTHS,
                                            n_sites=N_SITES, n_samples=N_SAMPLES,
                                            model=MODEL_PATH, rng_seed=RNG_SEED))
    m = load_numpy_model(MODEL_PATH)
    records, rep, ver = load_verified()
    gm = gene_map()
    symbols = sorted({v.upper() for v in gm.values()})
    wanted = {ALIASES.get(s, s) for s in symbols}
    print(f"Huesken genes: {len(symbols)}; querying GENCODE for {len(wanted)} symbols")
    tx = load_transcripts(GENCODE, wanted)
    print(f"GENCODE matched {len(tx)} gene symbols, {sum(len(v) for v in tx.values())} transcripts")

    # locate each siRNA's target site (the mRNA strand) inside its gene's transcripts
    sites, unmatched_gene, unmatched_site = [], 0, 0
    for r in records:
        anti = r["sequence"]
        g = gm.get(anti)
        if not g: continue
        key = ALIASES.get(g.upper(), g.upper())
        if key not in tx:
            unmatched_gene += 1; continue
        target = rc(anti)
        hit = None
        for ti, seq in enumerate(tx[key]):
            pos = seq.find(target)
            if pos >= 0:
                hit = (key, ti, pos, len(seq)); break
        if hit is None:
            unmatched_site += 1; continue
        sites.append({"antisense": anti, "efficacy": r["efficacy"], "gene": key,
                      "tx": hit[1], "pos": hit[2], "tx_len": hit[3]})
    print(f"target sites located: {len(sites)}  "
          f"(gene absent from GENCODE: {unmatched_gene}, site not found in transcript: {unmatched_site})")
    if not sites:
        record(run, "f9_target_context", {"status": "no_sites"}); raise SystemExit(3)

    rng = np.random.RandomState(RNG_SEED)
    chosen = [sites[i] for i in rng.choice(len(sites), size=min(N_SITES, len(sites)), replace=False)]
    RNA.init_rand(RNG_SEED)

    out = {}
    for L in LENGTHS:
        rows = []
        for st in chosen:
            seq = tx[st["gene"]][st["tx"]]
            centre = st["pos"] + 10
            lo = max(0, centre - L // 2)
            win = seq[lo:lo + L]
            if len(win) < L or set(win) - set("ACGU"): continue
            n = len(win)
            X = generic_features(win)
            backbone = [(i, i + 1) for i in range(n - 1)]
            lat = sound_lattice(win, floor=0.0, canonical_only=True)
            _, hi, passes = exact_interval(m, n, backbone, lat["optional"], X)
            y_mfe = m.forward(X, build_A(n, backbone + sorted(pairs_of(lat["mfe_structure"])), [], []))
            ens_max, sets = sample_max(win, m, X, backbone, n, N_SAMPLES)
            opt = set(lat["optional"])
            cov = float(np.mean([P <= opt for P in sets]))
            rows.append(dict(gene=st["gene"], k=lat["k"], worst_case=hi, y_mfe=y_mfe,
                             ens_max=ens_max, coverage=cov,
                             sound=bool(ens_max <= hi + 1e-9),
                             mfe_gap=(ens_max - y_mfe) / max(abs(y_mfe), 1e-12),
                             slack=hi / max(ens_max, 1e-12), passes=passes,
                             lattice_log10=lat["lattice_size_log10"]))
        med = lambda k: float(statistics.median([r[k] for r in rows]))
        out[L] = dict(n=len(rows), median_k=med("k"), max_k=max(r["k"] for r in rows),
                      median_lattice_log10=med("lattice_log10"),
                      median_worst_case=med("worst_case"), median_y_mfe=med("y_mfe"),
                      median_ens_max=med("ens_max"), median_slack=med("slack"),
                      median_mfe_gap=med("mfe_gap"),
                      coverage=float(np.mean([r["coverage"] for r in rows])),
                      sound_all=all(r["sound"] for r in rows), rows=rows)
        o = out[L]
        print(f"\nwindow {L} nt  (n={o['n']} sites, {N_SAMPLES} Boltzmann samples each)")
        print(f"  k (optional canonical pairs) : median {o['median_k']:.0f}, max {o['max_k']}  "
              f"-> lattice 10^{o['median_lattice_log10']:.0f}, certified in {rows[0]['passes']} passes")
        print(f"  MFE point estimate           : {o['median_y_mfe']:.1f}")
        print(f"  sampled ensemble max         : {o['median_ens_max']:.1f}")
        print(f"  SOUND certified worst case   : {o['median_worst_case']:.1f}")
        print(f"  ensemble coverage            : {o['coverage']*100:.1f}%  "
              f"sound on every site: {o['sound_all']}")
        print(f"  MFE understates ensemble by  : {o['median_mfe_gap']*100:.1f}%   "
              f"certificate slack vs sampled max: {o['median_slack']:.1f}x")

    record(run, "f9_target_context", {"by_length": out, "n_sites_located": len(sites),
                                      "genes_matched": sorted(tx), "model": MODEL_PATH})

if __name__ == "__main__": main()
