# Artifact reproduction audit

Audit date **2026-09-05**. Subject: the results narrative dated 4 September 2026
(`RESULTS.md`, and the shared report derived from it).

Scope: verify that every headline number is traceable to a script, a configuration, an
input, an output record, a code state and a software version; re-run the reproduction
commands that exist; and list every claim in the repository that is stale or withdrawn.

**No code was changed.** This audit adds documentation only. Experiment runs write to
`runs/`, which is gitignored, so the reproductions below leave the tracked tree clean.

---

## 1. Repository state at audit time

| item | value |
|---|---|
| HEAD | `35c0716` kappa C7-C9: the relaxation gap is a walk-count ratio |
| working tree | clean before and after this audit (tracked files) |
| tags | `prereg` → `2321b11` (scaffold + frozen pre-registration), the only tag |
| commits | 12, all dated 2026-09-04 |
| tracked files | 39 |
| run records | 38 pre-audit, plus 4 written by this audit |

### Software, identical across all 38 pre-audit runs and all audit runs

| component | version |
|---|---|
| Python | 3.13.12 |
| numpy | 2.4.2 |
| ViennaRNA (`RNA`) | 2.6.4 |
| torch | 2.11.0+cu130 |
| platform | Linux-6.17.0-22-generic-x86_64-with-glibc2.39 |

### Commit → experiment map

| commit | introduces |
|---|---|
| `2321b11` | scaffold, pre-registration, `f1`, `f2`, `f3`, `f4`, core `certmp` modules |
| `87a3aa9` | first results narrative for F1–F4 |
| `03e2bc6` | `f5_stress_extremality.py` |
| `92be4f8` | `f3b_relaxation_slack.py` |
| `72a0876` | `f6_trained.py`; last change to `certmp/data.py` |
| `45a72f0` | `make data` fetch target |
| `02c0b14` | `f7_soundness.py`, `certmp/train.py`, `experiments/_huesken.py`; last change to `certmp/certify.py`, `f3`, `f4` |
| `9101597` | `f8_baseline.py` |
| `b6f06d4` | `f9_target_context.py`; last change to `certmp/ensemble.py` and `f7_soundness.py` |
| `08428a5` | `f10_maximal.py`, `f10b_convergence.py`, `certmp/maximal.py` |
| `f842492` | `c3`–`c6`, `certmp/aggregators.py`; last change to `certmp/models.py`, `certmp/reach.py` |
| `35c0716` | `c7`–`c9`, `certmp/kappa.py` |

### Inputs (gitignored, third-party, not vendored)

| file | sha256 (first 16) | bytes |
|---|---|---|
| `data/TrainAll2182.txt` | `4919d99a6d5bffe2` | 61,096 |
| `data/TestAll249.txt` | `7ebfc59a0d0c2116` | 6,972 |
| `data/Huesken_2431_annotated.tsv` | `9e9be35db9c6a092` | 266,770 |
| `data/gencode/gencode.v47.pc_transcripts.fa.gz` | `5fed4bff68c46744` | 48,625,053 |
| `data/gencode_windows.txt` | `5a9cb5414c9e514d` | 3,636 |
| `data/models/monotone_generic.json` | `a696d123f396631e` | 31,512 |

`data/Huesken_2431_annotated.tsv`, `data/gencode_windows.txt` and
`data/models/monotone_generic.json` are **derived** artifacts with no fetch target in the
`Makefile` — see [§6 D-07](#6-discrepancy-table).

---

## 2. How provenance is recorded, and two defects in it

`certmp/provenance.py` writes, per run: experiment name, timestamp, `git_sha`,
`git_dirty`, Python version, platform, library versions, the config dict and a
`config_sha256`. That is more than most repositories carry. Two defects limit what it can
actually establish.

### D-01 · 33 of 38 pre-audit runs were made from a dirty tree

Only five pre-audit runs recorded `git_dirty: false`. For every other run, the recorded
`git_sha` is the commit that was checked out at the time — which, in this workflow, is the
commit *before* the experiment was committed. **The recorded SHA therefore does not
identify the code that ran.** Examples:

| run | experiment | recorded SHA | commit that actually contains the code |
|---|---|---|---|
| `20260904-020229` | `f9_target_context` | `9101597` (f8) | `b6f06d4` |
| `20260904-023117` | `f10_maximal` | `b6f06d4` (f9) | `08428a5` |
| `20260904-043511` | `c8_kappa_gencode` | `f842492` (C1–C6) | `35c0716` |

In each case the recorded SHA predates the existence of the script named in the same
record. The mapping in §3 gives the *correct* code state for each number; it was
reconstructed with `git log --diff-filter=A` and `git log -1` per file, not read from the
provenance records.

### D-02 · Two runs share a SHA and a config hash but produced different numbers

`runs/20260904-020117` and `runs/20260904-020253` both record `git_sha 9101597`,
`git_dirty true`, `config_sha256 52cff106bced46e7`. Their results differ materially:

| run | floor 0, canonical: bound ÷ ensemble max |
|---|---|
| `20260904-020117` | 50.00 |
| `20260904-020253` | **29.25** ← the value the report uses |

The difference is the backbone-as-mandatory fix, applied to `f7_soundness.py` between the
two runs while the tree was dirty. The recorded provenance cannot distinguish them. Two
contributing causes:

- `config_sha256` hashes only the config dict, so it is not unique to an experiment
  (`f1_extremality` and `c3_aggregator_characterization` share hash `f2750dc45dfe7bb0`
  because both configs are `{trials: 25, n: 7, k: 8}`).
- F7's config records `floors` as `[0.0, 0.0, 1e-06, …]` without the paired
  `canonical_only` flag, so the config alone does not describe the sweep.

**Consequence for the audit:** for every superseded run, authority comes from the
narrative in `RESULTS.md` plus the file's last-modifying commit, not from the run record.
Each row in §3 names the run that the report's number actually came from.

---

## 3. Headline number → source map

Column *code state* is the commit that contains the code as run, corrected for D-01.
All runs used the software versions in §1. All configs are reproduced verbatim from the
run's `provenance.json`.

### F1 · extremality, registered predictions R1–R6

| number | value | source |
|---|---|---|
| R1 sum, exact trials | 25/25, gap 0.000e+00 | |
| R2 max, exact trials | 25/25, gap 0.000e+00 | |
| R3 mean | 0/25, max gap 1.335e−01 | |
| R4 degree-norm | 0/25, max gap 8.853e−02 | |
| R5 signed weights | 11/25, max gap 1.108e+00 | |
| R6 negative features | 6/25, max gap 9.195e−01 | |

- script `experiments/f1_extremality.py` · command `make f1`
- config `{"trials": 25, "n": 7, "k": 8}` · hash `f2750dc45dfe7bb0`
- inputs none (synthetic, seeded `RandomState(0..24)`)
- output `runs/20260904-034428/f1_extremality.json`
- recorded SHA `08428a5` (dirty) · **code state `2321b11`** (file unchanged since scaffold)
- **reproduced 2026-09-05** → `runs/20260905-122427`, clean tree at `35c0716`, identical

### F2 · lattice size, registered prediction R7

| number | value |
|---|---|
| median k at 40 / 80 / 150 nt | 16.0 / 48.5 / 86.0 |
| median lattice at 150 nt | 10^25.9 |
| vacuous | false |

- script `experiments/f2_lattice_size.py` · command `make f2`
- config `{"lengths": [40,80,150], "n_seq": 30, "lo": 0.05, "hi": 0.9}` · hash `88235d01f24275bf`
- inputs none (random sequences, `random.seed(0)`); ViennaRNA 2.6.4 fold
- output `runs/20260904-020404/f2_lattice_size.json`
- recorded SHA `b6f06d4`, **`git_dirty: false`** · code state `2321b11`
- not re-run in this audit (not in the requested set)

### F3 · retired

Values (110.9 %, etc.) are retained in `RESULTS.md` under an explicit RETIRED heading.
Source `runs/20260904-015525/f3_certificate.json`. Superseded by F3b and F6. Not audited
as a live claim.

### F3b · relaxation slack

| number | value |
|---|---|
| MFE vs lattice worst case | 110.9 % |
| MFE vs sampled ensemble max | 13.4 % |
| share of gap that is slack | 87.9 % |
| median slack ratio | 1.82× |
| coverage, median / worst | 73.8 % / 43.1 % |
| lower-bound violations | 2 of 20 sequences |
| MFE pairs inside lattice | 20/20, 0 failures |

- script `experiments/f3b_relaxation_slack.py` · command `make f3b`
- config `{"n_samples": 1000, "n_seq": 20, "length": 60, "lo": 0.05, "hi": 0.9, "rng_seed": 20260904}` · hash `7cb64197412a3d43`
- inputs none (random sequences `random.seed(1)`, `RNA.init_rand(20260904)`)
- output `runs/20260904-034729/f3b_relaxation_slack.json`
- recorded SHA `08428a5` (dirty) · **code state `92be4f8`**
- earlier runs `…012711`, `…012737` used an unseeded RNG (config hash `4f9a47731d357904`) and are superseded

### F5 · stress audit (post-hoc)

| number | value |
|---|---|
| sum | 400 live / 400, 0 void, 0 violations |
| max | 286 live / 400, 114 void, 0 violations |
| max saturation, mandatory subgraph | live fraction 1.00 → 0.00 as density rises |
| RNA mandatory subgraph max degree | 1 (median density 0.28 %) |

- script `experiments/f5_stress_extremality.py` · command `make f5`
- config `{"trials": 400, "pre_registered": false}` · hash `efd8658e6b8825e5`
- inputs none for parts A/B; part C folds random 60 nt sequences
- output `runs/20260904-034701/f5_stress_extremality.json`
- recorded SHA `08428a5` (dirty) · **code state `03e2bc6`**
- **reproduced 2026-09-05** → `runs/20260905-122429`, clean tree at `35c0716`, identical

### F6 · trained model

| number | value |
|---|---|
| held-out Spearman / Pearson | 0.5916 / 0.5564 |
| gene-disjoint Spearman / Pearson | 0.5432 / 0.5307 (n = 290) |
| published baseline | Pearson 0.66 |
| gene leakage, published split | 30 of 30 genes shared |
| trained MFE vs lattice / vs ensemble | 12.9 % / 11.4 % |
| trained slack ratio | 1.00× (median k = 6) |
| torch↔numpy parity | 4.27e−09 |

- script `experiments/f6_trained.py` · command `make data && make f6`
- config `{"val_seed": 20260904, "epochs": 250, "lr": 0.01, "d_hid": 32, "n_layers": 2, "n_samples": 200, "batch": 128, "features": "one-hot (position, nucleotide)"}` · hash `5b0391880d68ccd7`
- inputs `data/TrainAll2182.txt`, `data/TestAll249.txt`, `data/Huesken_2431_annotated.tsv`
- output `runs/20260904-013609/f6_trained.json`
- recorded SHA `92be4f8` (dirty) · **code state `72a0876`**
- Huesken split verification **reproduced 2026-09-05**; full training not re-run (see §7)

### F7 · soundness and the coverage/tightness sweep

| floor | canonical | k | coverage | worst | bound ÷ ensemble max | violations |
|---|---|---|---|---|---|---|
| 0 | no | 1770 | 100.0 % | 100.0 % | 179.08 | 0 / 20000 |
| 0 | yes | 610 | 100.0 % | 100.0 % | **29.25** | 0 / 20000 |
| 1e−6 | yes | 414 | 100.0 % | 100.0 % | 15.14 | 0 / 20000 |
| 1e−4 | yes | 190 | 99.7 % | 99.3 % | 4.69 | 0 / 20000 |
| 1e−3 | yes | 120 | 98.2 % | 96.6 % | 2.57 | 0 / 20000 |
| 1e−2 | yes | 61 | 91.0 % | 79.2 % | 1.48 | 0 / 20000 |
| 0.05 | yes | 34 | 73.8 % | 43.1 % | 1.13 | 15 / 20000 |
| 0.10 | yes | 26 | 59.5 % | 36.3 % | 1.01 | 199 / 20000 |
| 0.20 | yes | 20 | 41.9 % | 8.1 % | 0.98 | 1017 / 20000 |

- script `experiments/f7_soundness.py` · command `make f7`
- config hash `52cff106bced46e7` · inputs `data/models/monotone_generic.json`
- output **`runs/20260904-020253/f7_soundness.json`** — the third of three F7 runs
- recorded SHA `9101597` (dirty) · **code state `b6f06d4`**
- superseded: `…015420` (no canonical column), `…020117` (pre-backbone-fix, 50.00 not 29.25). See D-02.
- figure `runs/20260904-020253/coverage_tightness.png`
- not re-run in this audit (not in the requested set)

### F8 · matched baseline

| encoding | model | Spearman | Pearson |
|---|---|---|---|
| positional | monotone | 0.592 ± 0.004 | 0.557 ± 0.004 |
| positional | unconstrained | 0.594 ± 0.013 | 0.588 ± 0.013 |
| generic | monotone | 0.389 ± 0.012 | 0.372 ± 0.010 |
| generic | unconstrained | 0.603 ± 0.007 | 0.594 ± 0.008 |
| positional, gene-disjoint | monotone | 0.544 ± 0.002 | 0.531 ± 0.001 |
| positional, gene-disjoint | unconstrained | 0.583 ± 0.023 | 0.570 ± 0.024 |

Cost of certifiability: +0.0025 (positional), +0.2141 (generic), +0.0388 (gene-disjoint).

- script `experiments/f8_baseline.py` · command `make f8`
- config `{"seeds": [0,1,2,3,4], "feature_sets": ["positional","generic"]}` · hash `cb552cb240c4f332`
- inputs `data/TrainAll2182.txt`, `data/TestAll249.txt`, `data/Huesken_2431_annotated.tsv`
- output `runs/20260904-015628/f8_baseline.json`
- recorded SHA `02c0b14` (dirty) · **code state `9101597`**
- not re-run in this audit (see §7)

### F9 · sound certificate on GENCODE target sites

| window | k | lattice | MFE | ensemble max | worst case | slack | coverage | MFE gap |
|---|---|---|---|---|---|---|---|---|
| 50 nt | 408 | 10^123 | 1.0 | 1.1 | 23.7 | 21.5× | 100 % | 4.4 % |
| 100 nt | 1744 | 10^525 | 1.5 | 1.5 | 146.3 | 94.3× | 100 % | 4.5 % |
| 150 nt | 3932 | 10^1183 | 1.9 | 2.0 | 446.4 | 222.8× | 100 % | 5.1 % |

Sites located 1720 of 2431; 28 of 30 gene symbols matched.

- script `experiments/f9_target_context.py` · command `make gencode && make f9`
- config `{"lengths": [50,100,150], "n_sites": 24, "n_samples": 300, "rng_seed": 20260904}` · hash `6b85eccfdd4bea68`
- inputs `data/gencode/gencode.v47.pc_transcripts.fa.gz`, `data/models/monotone_generic.json`, Huesken tables
- output `runs/20260904-020229/f9_target_context.json`
- recorded SHA `9101597` (dirty) · **code state `b6f06d4`**
- **reproduced 2026-09-05** → `runs/20260905-122512`, clean tree at `35c0716`, identical on all 24 reported values

### F10 · maximal-element route

| window | M | M95 | M99 | ratio 95 | ratio 99 | covered closure (95) | lattice slack |
|---|---|---|---|---|---|---|---|
| 50 nt | 361 | 79 | 264 | 1.008 | 1.020 | 10^4.82 – 10^6.49 | 21.3× over 10^123 |
| 100 nt | 1138 | 719 | 1046 | 1.027 | 1.031 | 10^10.24 – 10^12.48 | 91.0× over 10^521 |
| 150 nt | 2956 | 2518 | 2872 | 1.045 | 1.045 | 10^15.65 – 10^18.03 | 216.6× over 10^1154 |

Slack scaling fit: exponent 1.9284, intercept −2.2973, R² 0.9716, n = 36.

- script `experiments/f10_maximal.py` · command `make f10`
- config `{"lengths": [50,100,150], "n_samples": 5000, "n_windows": 12, "targets": [0.95,0.99], "rng_seed": 20260904}` · hash `00d716df56c5ae66`
- inputs `data/gencode/…fa.gz`, `data/models/monotone_generic.json`
- output `runs/20260904-023117/f10_maximal.json`
- recorded SHA `b6f06d4` (dirty) · **code state `08428a5`**

### F10b · convergence of the maximal set

M95 growth exponent in sample size: 50 nt `[0.226, 0.040, −0.011, 0.222]`; 150 nt
`[0.754, 0.855, 0.902, 0.751]`.

- script `experiments/f10b_convergence.py` · command `make f10b`
- config `{"lengths": [50,150], "sample_sizes": [1000,5000,20000], "n_windows": 4}` · hash `c3dd982f43396cc4`
- output `runs/20260904-023333/f10b_convergence.json`
- recorded SHA `b6f06d4` (dirty) · **code state `08428a5`**

### C3 · aggregator characterisation

sum / max / logsumexp isotone, max at full, 25/25 each. min antitone, max at **empty**,
25/25. mean / degnorm / std neither, 0/25 each. `characterization_holds: true`.

- script `experiments/c3_aggregator_characterization.py` · command `make c3`
- config `{"n": 7, "k": 8, "trials": 25}` · hash `f2750dc45dfe7bb0` (collides with F1, see D-02)
- output `runs/20260904-034519/c3_aggregator_characterization.json`
- recorded SHA `08428a5` (dirty) · **code state `f842492`**

### C4 · affinity

sum+relu additivity violation **2.157e−16** (affine); all 12 aggregation × activation
configurations agreed with prediction.

- script `experiments/c4_affinity.py` · command `make c4`
- config `{"trials": 30, "n": 6, "d_in": 3}` · hash `e990f660cd37a30a`
- output `runs/20260904-034529/c4_affinity.json`
- recorded SHA `08428a5` (dirty) · **code state `f842492`**

### C5 · slack scaling

Exponent per layer 0.912 at depths 1–4, R² 0.9996 throughout. Unbounded slack: 9.0 (m=3),
35.8, 143.2, 572.6, 2290.3 (m=48). Join-closed families exactly 1.0, brute-force verified.

- script `experiments/c5_slack_scaling.py` · command `make c5`
- config `{"m_range": [3,12], "depths": [1,2,3,4]}` · hash `87deef67304b9992`
- output `runs/20260904-034559/c5_slack_scaling.json`
- recorded SHA `08428a5` (dirty) · **code state `f842492`**

### C6 · closure

52/52 live configurations exact, 0 violations, 2 voids excluded, control discriminates
(4/4 broken by a non-monotone activation).

- script `experiments/c6_closure.py` · command `make c6`
- config `{"n": 6, "k": 7, "trials": 15, "d": 8}` · hash `f5a0727881502ce1`
- output `runs/20260904-034625/c6_closure.json`
- recorded SHA `08428a5` (dirty) · **code state `f842492`**
- **reproduced 2026-09-05** → `runs/20260905-122446`, clean tree at `35c0716`, identical

### C7 · kappa characterisation

K1–K4 all pass; max relative error **2.121e−16**. Chain (matching → deg≤2 → deg≤3 →
all-subsets): L=1 `2.3846 / 1.8235 / 1.4762 / 1.0`; L=2 `6.28 / 3.3404 / 2.0933 / 1.0`;
L=3 `15.5714 / 5.5693 / 2.6961 / 1.0`.

- script `experiments/c7_kappa_characterization.py` · command `make c7`
- config `{"n": 7, "n_candidate": 12, "seeds": 8, "depths": [1,2,3]}` · hash `e9637ab28a1e1146`
- output `runs/20260904-043204/c7_kappa_characterization.json`
- recorded SHA `f842492` (dirty) · **code state `35c0716`**

### C8 · kappa against measured slack

| window | measured | κ₁ | κ₂ | measured ÷ κ₂ | bounded |
|---|---|---|---|---|---|
| 50 nt | 21.5 | 10.69 | 118.5 | 0.181 | yes |
| 100 nt | 94.3 | 21.95 | 494.1 | 0.191 | yes |
| 150 nt | 222.8 | 32.37 | 1098.0 | 0.203 | yes |

Growth: measured ×4.39, ×2.36 vs κ₂ ×4.17, ×2.22. `real_windows: true`.

- script `experiments/c8_kappa_gencode.py` · command `make c8`
- config `{"windows": "data/gencode_windows.txt", "lengths": [50,100,150], "sample_sweep": [500,2000,5000], "depths": [1,2,3], "max_seqs": 12, "floor": 0.0}` · hash `00fbfa1ccacb8c60`
- inputs `data/gencode_windows.txt` (derived, see D-07)
- output `runs/20260904-043511/c8_kappa_gencode.json`
- recorded SHA `f842492` (dirty) · **code state `35c0716`**
- the `measured` column is **hard-coded** in the script as `MEASURED = {50: 21.5, 100: 94.3, 150: 222.8}`, not read from F9's record — see D-06

### C9 · boundary of the characterisation

Equality holds only for uniform features + sum + relu; `kappa_bounds_all_cases: true`
across varied features, max aggregation and tanh.

- script `experiments/c9_kappa_boundary.py` · command `make c9`
- config `{"n": 7, "n_candidate": 12, "seeds": 8}` · hash `8f48ec236e3da079`
- output `runs/20260904-043246/c9_kappa_boundary.json`
- recorded SHA `f842492` (dirty) · **code state `35c0716`**

---

## 4. Reproduction log, 2026-09-05

Only the smallest **official** command was used in each case; nothing was invented. All
four ran from a clean tree at `35c0716`, so their run records are unambiguous
(`git_dirty: false`).

| target | command | result | new run record |
|---|---|---|---|
| F1 | `make f1` | **identical** — all six rows match to the printed precision | `runs/20260905-122427` |
| F5 | `make f5` | **identical** — sum 400/400/0 violations; max 286 live, 114 void, 0 violations; part C max degree 1 | `runs/20260905-122429` |
| C6 | `make c6` | **identical** — 52/52, 0 violations, 2 voids, control 4/4 breaks | `runs/20260905-122446` |
| Huesken split | `make data` | **identical** — re-downloaded from `biodev.cea.fr/DSIR`, byte-identical sha256 to the pre-audit files; `{n_train: 2182, n_test: 249, overlap: 0}`, `verified: true`, `problems: []` | n/a (no run record written) |
| target-window certificate | `make f9` | **identical** — all 24 reported values across 50/100/150 nt match, including k, lattice exponent, worst case, slack and 100 % coverage | `runs/20260905-122512` |

Determinism note: F1, F5 and C6 are seeded numpy only. `make f9` involves ViennaRNA
partition-function folding and `pbacktrack` sampling seeded via `RNA.init_rand(20260904)`;
it reproduced exactly on the same ViennaRNA build. Reproduction across a different
ViennaRNA version is **not** established.

---

## 5. Stale and withdrawn claim scan

Searched the whole tree (excluding `runs/`) for: `top-k`, `topk`, `110.9`, `1.82`,
`safety`, `off-target`, `every structure`, `independent of k`, and old probability floors
(`lo=0.05`, `hi=0.90`, band language).

**`safety`** — no occurrences anywhere. Clean.

**`independent of k`** — 3 occurrences, all correct. `RESULTS.md:468` and
`PREREGISTRATION.md:80` state the C2-narrowed claim (two evaluations independent of k
against 2^k enumeration), which is accurate. `RESULTS.md:88` is an unrelated usage (the
max void rate is independent of k). No action.

**`110.9`** — 7 occurrences, all in `RESULTS.md`, all inside the retired-F3 section or
explicitly labelled as the superseded value. The retirement notice at `RESULTS.md:42`
precedes the first use. No action.

**`top-k` / `topk`** — correctly withdrawn in `certmp/certify.py` (raises
`NotImplementedError`), in `experiments/f4_topk_stability.py` (prints a retirement
notice) and in `RESULTS.md`. **One stale occurrence:** `README.md:18` still advertises
the withdrawn certificate. See D-03.

**`1.82`** — occurrences in `RESULTS.md` are correct in context (F3b's banded-lattice
figure, explicitly superseded) and the `1.824`/`1.823` hits are the unrelated C5 exponent
and C7 chain value. **One stale occurrence:** `README.md:11` presents 1.82× as the
current looseness. See D-04.

**`off-target`** — 1 occurrence, `README.md:9`. Substantive problem: see D-05.

**`every structure`** — 9 occurrences. Eight are correct in context (the backbone is in
every structure; the floor-0 lattice contains every structure). **One stale:**
`certmp/ensemble.py:5` still says pairs above the band are "present in essentially every
structure", which is the assumption F3b falsified. See D-08.

**Old probability floors** — `bpp_lattice(seq, lo=0.05, hi=0.90, floor=1e-3)` remains the
default in `certmp/ensemble.py:17`, and is still used by `f2`, `f3`, `f3b` and `f6`. This
is **correct and intentional**: the module docstring now carries a "BANDED LATTICE, NOT
SOUND" warning and directs callers to `sound_lattice()`. No action.

---

## 6. Discrepancy table

Severity: **A** = affects a reported number or a soundness claim · **B** = affects
traceability · **C** = documentation drift.

| ID | Sev | Location | Discrepancy | Evidence |
|---|---|---|---|---|
| D-01 | B | all runs | 33 of 38 pre-audit runs have `git_dirty: true`; the recorded SHA predates the script it names, so it does not identify the code that ran | §2 |
| D-02 | B | `runs/20260904-020117`, `…020253` | Same `git_sha` **and** same `config_sha256`, materially different results (50.00 vs 29.25). `config_sha256` is not unique per experiment; F7's config omits `canonical_only` | §2 |
| D-03 | C | `README.md:18` | Advertises "threshold certificate + exact top-k rank stability". `certify_topk` was withdrawn and raises `NotImplementedError` | `certmp/certify.py:41-49` |
| D-04 | A | `README.md:11` | States the bound is "loose by a measured median factor of 1.82". That is the superseded **banded, unsound** lattice figure from F3b with an untrained model. The sound figure is **29.25×** at 60 nt (F7) and 222.8× at 150 nt (F9) | `runs/20260904-020253`, `runs/20260904-020229` |
| D-05 | A | `README.md:9` | Frames the certified quantity as "predicted **off-target** risk". The trained model is fit to Huesken **efficacy** — normalised inhibition of the *intended* target (`%Inhibition` in the annotated table). No off-target quantity is measured anywhere in the repository | `certmp/data.py`, `data/Huesken_2431_annotated.tsv` header |
| D-06 | B | `experiments/c8_kappa_gencode.py` | The `MEASURED = {50: 21.5, 100: 94.3, 150: 222.8}` comparison baseline is hard-coded, not read from F9's run record. Correct today, but it will silently disagree if F9 is re-run with different settings | script source |
| D-07 | B | `data/` | Three inputs are **derived with no fetch or build target**: `Huesken_2431_annotated.tsv` (produced by a join against siRNAEfficacyDB), `gencode_windows.txt` (extracted from GENCODE by an ad-hoc snippet), `models/monotone_generic.json` (trained by an ad-hoc snippet). `make data` and `make gencode` do not create them, so F7/F9/F10/C8 cannot be reproduced from a clean checkout | `Makefile`, §1 |
| D-08 | C | `certmp/ensemble.py:5` | Docstring still says pairs above the band are "present in essentially every structure". F3b measured a median 96.6 % (worst 83.0 %) of samples containing *all* mandatory pairs — the assumption that made the lower bound unsound | `runs/20260904-034729` |
| D-09 | A | `RESULTS.md:398` | F10 100 nt covered closure reported as "10^10 to 10^13". The run record gives 10^10.2350 – **10^12.4784**, which rounds to 10^12, not 10^13. Propagated into the shared report | `runs/20260904-023117/f10_maximal.json` |
| D-10 | C | `README.md:14-20` | Module listing omits `aggregators.py`, `train.py`, `maximal.py`, `kappa.py`, and describes a two-sided certificate. Predates commits `02c0b14` onward | `git ls-files certmp/` |
| D-11 | C | `README.md:2-4` | "provided the aggregation is sum or max" — superseded by C3, which admits any aggregator monotone under multiset inclusion, including `logsumexp` and `min` | `runs/20260904-034519` |
| D-12 | B | F9 vs F10 lattice sizes | F9 reports 10^525 / 10^1183 at 100 / 150 nt, F10 reports 10^521 / 10^1154 for the same lengths. **Not an error** — F9 uses 24 sites, F10 uses 12 — but the report places the two side by side without stating that they are different window sets | both run configs |

**Nothing in `RESULTS.md` was edited, and no artifact was regenerated.** D-09 and the
README items are recorded here for the maintainer to action deliberately.

---

## 7. Unreproduced items

Marked unreproduced rather than assumed. None was skipped for lack of a command; the
commands exist and the blockers are stated.

| item | status | blocker / reason |
|---|---|---|
| F2, F3b, F7, F8, F10, F10b, C3, C4, C5, C7, C8, C9 | **not re-run** | Outside the set requested for this audit. All have official `make` targets and recorded outputs; each is traceable per §3 |
| F6 full training | **partially reproduced** | The dataset half (`make data`) reproduced byte-identically and re-verified. Training was not re-run: F6 needs `data/Huesken_2431_annotated.tsv`, which has no fetch target (D-07) |
| Any run from a clean checkout | **blocked** | D-07. `make data` and `make gencode` fetch two of six inputs. The other three are derived by snippets that exist only in session history, not in the repository |
| Cross-version reproduction | **not attempted** | All 42 runs used one ViennaRNA build (2.6.4) and one numpy (2.4.2). Sensitivity of the sampled quantities to a different ViennaRNA is unknown |
| `runs/20260904-020253/coverage_tightness.png` | **not regenerated** | Regenerating would require re-running F7; the source record is intact |

---

## 8. Acceptance criteria

| criterion | status |
|---|---|
| Every reported result reproduced with provenance, or marked unreproduced | **met** — §3 traces all 17 experiments; §4 reproduces the 5 requested; §7 marks the rest |
| No code changes | **met** |
| `git diff` contains documentation only | **met** — the only change is this file, `docs/artifact_reproduction.md` |

Verify with:

```
git status --porcelain      # docs/artifact_reproduction.md only
git diff --stat HEAD -- certmp experiments tests Makefile   # empty
```
