#!/usr/bin/env python3
"""Derived statistics for the lipid restoration runs, from their per-seed arrays.

Reads the two results files written by lipid_restoration.py and writes
derived_statistics.json. Adds nothing the runs did not record: every field here is a
function of the per-seed held-out contrast MSEs, the recorded contrast energy and the
two floors.

    python3 derive_statistics.py > derived_statistics.json
"""
import json, sys
import numpy as np
from math import comb

FILES = (("seeds20", "lipid_restoration_results.json"),
         ("seeds40", "lipid_restoration_results_40.json"))

def sign_test(diff):
    """Exact two-sided binomial test on the per-seed sign."""
    n = len(diff); k = int((diff > 0).sum())
    return k, sum(comb(n, i) for i in range(n + 1)
                  if min(i, n - i) <= min(k, n - k)) / 2 ** n

out = {}
for tag, f in FILES:
    d = json.load(open(f))
    r2 = np.array(d["heldout_r2"]); r8 = np.array(d["heldout_r8"])
    diff = r2 - r8; E = d["energy"]
    k, p = sign_test(diff)
    out[tag] = dict(
        seeds=len(diff), capped_fits=d["capped_fits"], max_iter=d["max_iter"],
        folds=d["folds"], q=d["q"],
        mean_heldout_contrast_mse_radius2=float(r2.mean()),
        mean_heldout_contrast_mse_radius8=float(r8.mean()),
        paired_mean_radius2_minus_radius8=float(diff.mean()),
        bootstrap_ci_95=d["ci"],
        seeds_favouring_radius8=k, sign_test_two_sided_p=p,
        trimmed_mean_drop_two_largest_abs=float(np.sort(diff)[1:-1].mean()),
        zero_interaction_baseline=E,
        share_of_energy_recovered_radius2=float((E - r2.mean()) / E),
        share_of_energy_recovered_radius8=float((E - r8.mean()) / E),
        seeds_worse_than_predicting_zero_radius2=int((r2 > E).sum()),
        seeds_worse_than_predicting_zero_radius8=int((r8 > E).sum()),
        gap_over_radius2_floor=float(diff.mean() / d["floor_r2"]),
        floor_radius2=d["floor_r2"], floor_radius8=d["floor_r8"],
        rank_CH_radius2=d["rank_r2"], classes_radius2=d["classes_r2"])

out["definitions"] = dict(
    contrast_mse="q^-1 ||I - C p||^2 over the q=800 anchored tail-by-tail rectangles, "
                 "p the held-out prediction vector",
    paired_difference="radius-2 contrast MSE minus radius-8 contrast MSE, per seed; "
                      "positive means the lossless encoder had the lower error",
    zero_interaction_baseline="q^-1 ||I||^2, the contrast MSE of predicting no interaction",
    bootstrap_ci_95="as written by lipid_restoration.py: 20,000 resamples of the per-seed "
                    "differences at its own fixed seed, a percentile interval on the mean",
    sign_test="exact two-sided binomial test on the per-seed sign, not a test on the mean",
    derived_from="the per-seed arrays in the two results files in this directory")

json.dump(out, sys.stdout, indent=1)
sys.stdout.write("\n")
