# Restoring the lossy lipid encoder: held-out interaction error

This directory holds one experiment and nothing from the paper. It is not referenced by the
manuscript, no reported value depends on it, and it was run after the revision was built.

## What was compared

The AGILE lipid panel of Appendix A10.23, 1,100 single-compound cells over one head-group factor
and two tail factors, with the 800 anchored tail-by-tail interaction contrasts of that appendix.
Two encoders, same panel, same rows, same folds:

- **radius 2** Morgan count fingerprints, 2,048 bins: merges the 1,100 states into 920 classes,
  rank(CH) 740 against rank(C) 800, encoding floor 0.3150033148.
- **radius 8** Morgan count fingerprints, 2,048 bins: separates all 1,100 states, rank(CH) 800,
  floor exactly zero.

A 64-unit MLP is fitted on each arm and scored on held-out interaction contrasts, 5-fold grouped by
head group, so each fold predicts the tail-by-tail interactions of head groups it never saw.
Training rows are identical across the two arms, so training support is held fixed and only
representability differs. The reported quantity is contrast MSE, `q^-1 ||I - C p||^2` over the 800
rectangles.

## What came out

Positive means the lossless encoder had the lower error.

| run | paired mean | bootstrap 95% CI | favours radius 8 | sign-test p | fits at the iteration cap |
|---|---|---|---|---|---|
| 20 seeds | +0.536458 | [+0.050279, +1.059340] | 13 of 20 | 0.263 | 0 of 200 |
| 40 seeds | +0.658342 | [+0.326735, +0.997114] | 28 of 40 | 0.017 | 0 of 400 |

Restoring the encoding lowered held-out interaction error on this panel. Three qualifications
belong with that sentence, and all three are measured:

1. **Both arms sit near the no-interaction baseline.** Predicting no interaction at all gives
   contrast MSE 12.425324, the recorded contrast energy. Over 40 seeds radius 2 averages 12.207818
   and radius 8 11.549476, so they recover 1.75% and 7.05% of that energy. Radius 2 is worse than
   predicting zero in 16 of the 40 seeds, radius 8 in 11.
2. **The gap is larger than the floor explains.** The radius-2 floor is 0.315003 and the measured
   gap is 0.658342, 2.1 times it. The floor bounds what the merged classes permit; it did not
   predict the size of the improvement, in either direction.
3. **It is one panel, one architecture, one evaluation split, with no fitted model in the paper.**
   Nothing here transfers to modified siRNA, and nothing here is a claim about any published
   pipeline run end to end.

## Reproducing it

The panel file is already in the package at
`results/analyses/independent_lossy_encoder/inputs/lantern.csv`, sha256
`5db98d48f635f939196c9125b5ee15f217920b75cef2f8ccbe15e1c1dbb13d45`, which is the `panel_sha256`
recorded in `diagnostic_results.json`.

```
python3 lipid_restoration.py <path to>/inputs/lantern.csv --seeds 40
python3 derive_statistics.py > derived_statistics.json
```

`lipid_restoration.py` refuses to fit unless its self-check reproduces the published diagnostic
first: 920 classes, rank(CH) 740, floor 0.3150033148 at 2.5352% of energy, contrast energy
12.425324. Both runs here passed that check, and both logs begin with it.

Needs rdkit, numpy and scikit-learn. These runs used **rdkit 2026.03.1**, numpy 2.4.2,
scikit-learn 1.8.0, scipy 1.18.0 on CPython 3.13.12. The rdkit version is the one behind the
recorded diagnostic. Note that `lipid_restoration.py` builds its fingerprints without
`includeChirality`, which the recorded diagnostic passes; the self-check reproduces the floor to
ten digits either way, which is consistent with the partition being unchanged with and without
chirality as Appendix A10.23 states.

## Files

| file | what it is |
|---|---|
| `lipid_restoration.py` | the experiment, as run, including the self-check it aborts on |
| `run_seeds20.log.txt` | full stdout of the 20-seed run, per seed |
| `run_seeds40.log.txt` | full stdout of the 40-seed run, per seed |
| `lipid_restoration_results.json` | the 20-seed record, every per-seed value |
| `lipid_restoration_results_40.json` | the 40-seed record, every per-seed value |
| `derive_statistics.py` | computes the derived numbers above from those two records |
| `derived_statistics.json` | its output: means, intervals, sign test, baseline shares, floor ratio |

An earlier prototype of this experiment, run elsewhere with fits that hit the iteration cap, gave
the opposite sign: radius 2 lower by 0.5878 over 10 seeds, favouring radius 8 in 3 of 10. On the
same 10 seeds with nothing capped the difference is +0.359451, favouring radius 8 in 6 of 10. That
earlier run is not reproduced here and its records are not in this package.
