# Reproducible environment manifest

Captured **2026-09-05** from the interpreter that produced every run record under `runs/`.
Companion lock file: `docs/requirements-lock.txt` (226 packages, sha256 `c0132b85389ea9d0`).

---

## 1. Platform

| item | value |
|---|---|
| OS | Ubuntu 24.04.4 LTS (Noble Numbat) |
| kernel | Linux 6.17.0-22-generic, x86_64 |
| Python | 3.13.12, packaged by conda-forge, built with GCC 14.3.0 |
| accelerator | 2 CUDA devices visible; `torch.cuda.is_available() == True` |
| CUDA usage | **none of the logged results used the GPU.** `certmp/train.py` builds CPU tensors and never calls `.cuda()`; the certificate path is numpy-only |

---

## 2. Package versions, and a discrepancy that matters

| package | `__version__` (what provenance logged) | distribution version (pip) | used by |
|---|---|---|---|
| numpy | 2.4.2 | 2.4.2 | everything |
| **ViennaRNA / `RNA`** | **2.6.4** | **2.7.2** | `ensemble.py`, F2/F3b/F7/F9/F10/C8 |
| torch | 2.11.0+cu130 | 2.11.0 | `train.py`, F6/F8 only |
| matplotlib | 3.10.8 | 3.10.8 | F7 figure only |
| scipy | 1.18.0 | 1.18.0 | **installed but never imported** — `data.py` implements `spearman`/`pearson` locally |

### E-1 · The ViennaRNA version in every run record is not the installed version

`certmp/provenance.py` records `__import__(mod).__version__`. For `RNA` that returns
**2.6.4**, while `importlib.metadata.version("ViennaRNA")` returns **2.7.2**. All 42 run
records therefore log `"RNA": "2.6.4"`, which does not identify the installed
distribution — the attribute appears to report a bundled library version rather than the
package release.

**Why this is not cosmetic.** F7, F9, F10, F10b and C8 depend on ViennaRNA's partition
function and on `pbacktrack` sampling. Someone reconstructing the environment from a run
record would `pip install ViennaRNA==2.6.4` and obtain a **different** package from the
one that produced these numbers. Sampling behaviour across that boundary is untested.

**Action required (code change, not made here):** record
`importlib.metadata.version(dist)` alongside `__version__`, and re-record the mapping for
existing runs. Until then, treat `RNA 2.6.4` in any run record as meaning *"the RNA module
whose `__version__` string is 2.6.4, shipped in ViennaRNA distribution 2.7.2"*.

### E-2 · `requirements.txt` is unpinned and incomplete

```
numpy>=1.24
ViennaRNA>=2.6
torch>=2.0
```

Three problems: the bounds are open, so a fresh install does not reproduce the logged
environment; `matplotlib` is imported by `experiments/f7_soundness.py` but is **not
declared**; and `ViennaRNA>=2.6` would accept 2.7.2 or 2.6.4 indifferently, which E-1
shows are distinguishable. `docs/requirements-lock.txt` is the reproducible artifact;
`requirements.txt` is not.

---

## 3. Source state

Repository HEAD `35c0716`, tracked tree clean. sha256 (first 12) of every tracked source
file at capture time:

| file | sha256 |
|---|---|
| `Makefile` | `2cf3c0c07def` |
| `requirements.txt` | `2ce2b4adaa0b` |
| `certmp/__init__.py` | `ac84119dbb4b` |
| `certmp/aggregators.py` | `f09e2c0ee50a` |
| `certmp/certify.py` | `76bf4d979a8e` |
| `certmp/data.py` | `c818e718b8e9` |
| `certmp/ensemble.py` | `dd880623df6e` |
| `certmp/kappa.py` | `658253ebeb07` |
| `certmp/maximal.py` | `486029d141b1` |
| `certmp/models.py` | `75e82922ce39` |
| `certmp/provenance.py` | `bf89113cc54d` |
| `certmp/reach.py` | `0032bc35489f` |
| `certmp/train.py` | `e27789842ca1` |
| `certmp/void.py` | `7fb34ca0a46f` |
| `experiments/_huesken.py` | `b6ab3819335d` |
| `experiments/c3_aggregator_characterization.py` | `2aa715d00f50` |
| `experiments/c4_affinity.py` | `7b6999671d15` |
| `experiments/c5_slack_scaling.py` | `1aec2004268d` |
| `experiments/c6_closure.py` | `921285fb46c8` |
| `experiments/c7_kappa_characterization.py` | `51bcb8243169` |
| `experiments/c8_kappa_gencode.py` | `a536563167f0` |
| `experiments/c9_kappa_boundary.py` | `0fe78ef7d1aa` |
| `experiments/f1_extremality.py` | `9d882b3d3af9` |
| `experiments/f2_lattice_size.py` | `23783c200732` |
| `experiments/f3_certificate.py` | `bcd12a5351ec` |
| `experiments/f3b_relaxation_slack.py` | `ddd10574cbb9` |
| `experiments/f4_topk_stability.py` | `e47c69714253` |
| `experiments/f5_stress_extremality.py` | `0355457cf9ac` |
| `experiments/f6_trained.py` | `5dc075bfe27c` |
| `experiments/f7_soundness.py` | `1d88d0bcc076` |
| `experiments/f8_baseline.py` | `d02a079554b3` |
| `experiments/f9_target_context.py` | `36bec1c79ced` |
| `experiments/f10_maximal.py` | `233c2ccda3f6` |
| `experiments/f10b_convergence.py` | `81278860f50e` |
| `tests/test_sanity.py` | `b8fa993c7804` |

---

## 4. Data inputs

Third-party, gitignored, not vendored.

| file | sha256 (first 16) | bytes | how obtained |
|---|---|---|---|
| `data/TrainAll2182.txt` | `4919d99a6d5bffe2` | 61,096 | `make data` → biodev.cea.fr/DSIR |
| `data/TestAll249.txt` | `7ebfc59a0d0c2116` | 6,972 | `make data` → biodev.cea.fr/DSIR |
| `data/gencode/gencode.v47.pc_transcripts.fa.gz` | `5fed4bff68c46744` | 48,625,053 | `make gencode` → ftp.ebi.ac.uk |
| `data/Huesken_2431_annotated.tsv` | `9e9be35db9c6a092` | 266,770 | **no build target** |
| `data/gencode_windows.txt` | `5a9cb5414c9e514d` | 3,636 | **no build target** |
| `data/models/monotone_generic.json` | `a696d123f396631e` | 31,512 | **no build target** |

### E-3 · Three inputs have no build target

`make data` and `make gencode` fetch two of six inputs. The other three were produced by
ad-hoc snippets that exist only in session history. F7, F9, F10, F10b and C8 all depend on
at least one of them, so **those five experiments cannot be reproduced from a clean
checkout by anyone else.** This is the single largest reproducibility defect in the
project.

---

## 5. Reconstruction procedure

```bash
# 1. interpreter
conda create -n certmp python=3.13.12 && conda activate certmp

# 2. exact package set (NOT requirements.txt — see E-2)
pip install -r docs/requirements-lock.txt

# 3. source
git checkout 35c0716

# 4. fetchable inputs
make data      # verifies against published statistics, exits non-zero on mismatch
make gencode

# 5. verify the install reproduces the logged environment
python3 -c "import numpy,RNA,torch,importlib.metadata as m; \
  print(numpy.__version__, RNA.__version__, m.version('ViennaRNA'), torch.__version__)"
# expect: 2.4.2 2.6.4 2.7.2 2.11.0+cu130

# 6. smallest checks
make test      # 7 tests, ~1.2 s
make f1        # ~15 s
make c6        # ~20 s
```

Steps 1–6 reproduce F1, F2, F3b, F5, C3–C7, C9 and the Huesken verification. They do
**not** reproduce F6, F7, F8, F9, F10, F10b or C8, which need the three inputs in E-3.

---

## 6. Determinism

| source of randomness | control | verified |
|---|---|---|
| numpy | explicit `RandomState(seed)` per trial | yes — F1, F5, C6 reproduced bit-identically 2026-09-05 |
| Python `random` | `random.seed(0/1/2)` in F2, F3b, F3, F7, F9 | yes |
| ViennaRNA sampling | `RNA.init_rand(20260904)` in F3b, F7, F9, F10, C8 | yes — F9 reproduced identically on the same build |
| torch | `torch.manual_seed(seed)`, `Generator().manual_seed` in `train.fit` | not independently re-verified in this audit |
| GPU nondeterminism | not applicable — CPU only | — |
| cross-version determinism | **not established.** All runs used one ViennaRNA build | no |

---

## 7. Manifest defects, ranked

| id | defect | severity |
|---|---|---|
| **E-3** | three data inputs have no build target; five experiments unreproducible from a clean checkout | **high** |
| **E-1** | `RNA.__version__` (2.6.4) ≠ distribution version (2.7.2); every run record carries the ambiguous one | **high** |
| **E-2** | `requirements.txt` unpinned, and omits `matplotlib` | medium |
| — | `scipy` installed but unused; harmless, but it is in the lock and implies a dependency that does not exist | low |
