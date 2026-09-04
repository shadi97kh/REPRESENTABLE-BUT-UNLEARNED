"""F3b: quantify the relaxation slack that F3 leaves unmeasured.

F3 reports a gap between an MFE-only prediction and the lattice worst case. That gap
mixes two things: genuine Boltzmann-ensemble uncertainty, and slack from the fact that
the edge lattice is an OUTER relaxation. Most edge subsets are not valid nested
secondary structures, so the lattice maximum may be attained at a subset no real
structure realises. This experiment separates the two by sampling actual structures.

Method: RNA.pbacktrack draws structures from the Boltzmann distribution. Evaluate the
same model on each sampled structure and compare the true ensemble maximum against the
lattice maximum. The ratio is the slack.

Two soundness questions the pre-registration asserted rather than tested:
  Q1 Is every sampled structure inside the lattice, i.e. are its pairs a subset of
     mandatory union optional? Pairs below the 0.05 band are DISCARDED by bpp_lattice,
     so a structure using one is not covered and the upper bound does not apply to it.
  Q2 Does every sampled structure contain all mandatory pairs? The lower bound forces
     mandatory present. A structure missing one can score below the certified minimum,
     which would make the LOWER bound unsound.
Also verify the MFE structure itself lies inside the lattice.
"""
import random, statistics, numpy as np
import RNA
from certmp.models import MonoMPNN
from certmp.ensemble import bpp_lattice, onehot_features
from certmp.reach import exact_interval, build_A
from certmp.provenance import new_run, record

N_SAMPLES, N_SEQ, LENGTH = 1000, 20, 60
RNG_SEED = 20260904          # ViennaRNA pbacktrack is stochastic; seed it or results drift

def pairs_of(structure):
    """Dot-bracket -> set of 0-indexed base pairs."""
    stack, pairs = [], set()
    for i, c in enumerate(structure):
        if c == "(": stack.append(i)
        elif c == ")": pairs.add((stack.pop(), i))
    return pairs

def sample_structures(seq, n):
    md = RNA.md(); md.uniq_ML = 1          # required for vrna_pbacktrack
    fc = RNA.fold_compound(seq, md)
    _, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    return list(fc.pbacktrack(n))

def main():
    run = new_run("f3b_relaxation_slack",
                  dict(n_samples=N_SAMPLES, n_seq=N_SEQ, length=LENGTH, lo=0.05, hi=0.90, rng_seed=RNG_SEED))
    random.seed(1)                          # same sequences as F3
    RNA.init_rand(RNG_SEED)                 # reproducible Boltzmann sampling
    rows = []
    for t in range(N_SEQ):
        seq = "".join(random.choice("ACGU") for _ in range(LENGTH))
        lat = bpp_lattice(seq)
        if lat["k"] == 0: continue
        X = onehot_features(seq)
        m = MonoMPNN(X.shape[1], 8, 2, agg="sum", nonneg=True, seed=t)   # same model as F3
        n = lat["n"]
        in_lattice = set(lat["mandatory"]) | set(lat["optional"])
        mand = set(lat["mandatory"])

        lat_lo, lat_hi, _ = exact_interval(m, n, lat["mandatory"], lat["optional"], X)

        structs = sample_structures(seq, N_SAMPLES)
        ys, covered_ys, n_cov, n_has_mand = [], [], 0, 0
        for s in structs:
            P = pairs_of(s)
            y = m.forward(X, build_A(n, list(P), [], []))
            ys.append(y)
            cov = P <= in_lattice
            n_cov += cov
            n_has_mand += mand <= P
            if cov: covered_ys.append(y)

        mfe_pairs = pairs_of(lat["mfe_structure"])
        y_mfe = m.forward(X, build_A(n, list(mfe_pairs), [], []))
        row = dict(
            seq_id=t, k=lat["k"], n_samples=len(structs),
            n_unique=len(set(structs)),
            coverage=n_cov / len(structs),
            mandatory_present_frac=n_has_mand / len(structs),
            mfe_in_lattice=bool(mfe_pairs <= in_lattice),
            n_mandatory=len(mand),
            lattice_hi=lat_hi, lattice_lo=lat_lo,
            ens_max=max(ys), ens_min=min(ys),
            ens_max_covered=max(covered_ys) if covered_ys else None,
            slack_ratio=lat_hi / max(max(ys), 1e-12),
            upper_sound_covered=(max(covered_ys) <= lat_hi + 1e-9) if covered_ys else None,
            upper_sound_all=bool(max(ys) <= lat_hi + 1e-9),
            lower_sound_all=bool(min(ys) >= lat_lo - 1e-9),
            y_mfe=y_mfe,
            # what F3 reported: MFE vs the LATTICE worst case (contains slack)
            underestimate_lattice=(lat_hi - y_mfe) / max(abs(y_mfe), 1e-12),
            # the honest number: MFE vs the TRUE sampled-ensemble worst case
            underestimate_ensemble=(max(ys) - y_mfe) / max(abs(y_mfe), 1e-12),
        )
        rows.append(row)
        print(f"seq{t:02d} k={row['k']:3d} uniq={row['n_unique']:4d}  "
              f"cover={row['coverage']*100:5.1f}%  mand-in-struct={row['mandatory_present_frac']*100:5.1f}%  "
              f"ens_max={row['ens_max']:9.1f}  lat_hi={lat_hi:9.1f}  "
              f"slack x{row['slack_ratio']:5.2f}  MFE in lattice={row['mfe_in_lattice']}")

    med = lambda key: float(statistics.median([r[key] for r in rows]))
    summary = dict(
        n_seq=len(rows),
        median_slack_ratio=med("slack_ratio"),
        min_slack_ratio=min(r["slack_ratio"] for r in rows),
        max_slack_ratio=max(r["slack_ratio"] for r in rows),
        median_underestimate_lattice=med("underestimate_lattice"),
        median_underestimate_ensemble=med("underestimate_ensemble"),
        median_coverage=med("coverage"),
        min_coverage=min(r["coverage"] for r in rows),
        median_mandatory_present=med("mandatory_present_frac"),
        min_mandatory_present=min(r["mandatory_present_frac"] for r in rows),
        mfe_in_lattice_all=all(r["mfe_in_lattice"] for r in rows),
        mfe_in_lattice_failures=sum(not r["mfe_in_lattice"] for r in rows),
        upper_sound_covered_all=all(r["upper_sound_covered"] for r in rows if r["upper_sound_covered"] is not None),
        upper_sound_all_samples=all(r["upper_sound_all"] for r in rows),
        upper_sound_all_failures=sum(not r["upper_sound_all"] for r in rows),
        lower_sound_all_samples=all(r["lower_sound_all"] for r in rows),
        lower_sound_all_failures=sum(not r["lower_sound_all"] for r in rows),
    )
    ul, ue = summary["median_underestimate_lattice"], summary["median_underestimate_ensemble"]
    summary["slack_share_of_reported_gap"] = (ul - ue) / ul if ul else float("nan")
    print("\n--- summary ---")
    print(f"F3 reported gap, MFE vs LATTICE worst case   : {ul*100:.1f}%  (contains slack)")
    print(f"honest gap,      MFE vs SAMPLED ensemble max : {ue*100:.1f}%")
    print(f"share of F3's reported gap that is pure slack: {(ul-ue)/ul*100:.1f}%")
    print(f"median slack ratio lattice_hi / ensemble_max : x{summary['median_slack_ratio']:.2f} "
          f"(range x{summary['min_slack_ratio']:.2f} - x{summary['max_slack_ratio']:.2f})")
    print(f"Q1 sampled structures inside the lattice     : median {summary['median_coverage']*100:.1f}%, "
          f"worst {summary['min_coverage']*100:.1f}%")
    print(f"Q2 sampled structures containing ALL mandatory: median {summary['median_mandatory_present']*100:.1f}%, "
          f"worst {summary['min_mandatory_present']*100:.1f}%")
    print(f"MFE pairs inside mandatory u optional        : {len(rows)-summary['mfe_in_lattice_failures']}/{len(rows)} "
          f"({summary['mfe_in_lattice_failures']} failures)")
    print(f"upper bound sound on COVERED samples         : {summary['upper_sound_covered_all']}")
    print(f"upper bound sound on ALL samples             : {summary['upper_sound_all_samples']} "
          f"({summary['upper_sound_all_failures']} sequences violate)")
    print(f"lower bound sound on ALL samples             : {summary['lower_sound_all_samples']} "
          f"({summary['lower_sound_all_failures']} sequences violate)")
    record(run, "f3b_relaxation_slack", {"rows": rows, "summary": summary})

if __name__ == "__main__": main()
