# Claim ledger

Audit date **2026-09-05**. Every claim currently made anywhere in the repository or in the
4 September report, mapped to the function that computes it and the logged result that
supports it.

**Status vocabulary**

- **IMPLEMENTED** — code exists, ran, and a run record supports the stated number.
- **PROPOSED** — described in `docs/stage1_novelty_audit.md` as future work. No code.
- **UNSUPPORTED** — asserted somewhere in the repo but contradicted, superseded, or
  never measured.

Run records are under `runs/`; all used numpy 2.4.2, ViennaRNA 2.6.4, torch 2.11.0+cu130.

---

## A. Implemented and supported

| # | Claim | Function | Logged result | Status |
|---|---|---|---|---|
| A1 | Endpoint exactness for sum + non-negative weights + non-negative features | `reach.endpoint_gap` → `reach._endpoints` | `runs/20260904-034428` R1 25/25, gap 0.000e+00. Re-run `runs/20260905-122427` identical | IMPLEMENTED |
| A2 | Same for max | as A1 | same run, R2 25/25, gap 0.000e+00 | IMPLEMENTED |
| A3 | mean and degree-norm are **not** endpoint-exact | as A1 | R3 0/25 max gap 1.335e−01; R4 0/25 max gap 8.853e−02 | IMPLEMENTED |
| A4 | Signed weights break exactness | as A1 | R5 11/25 exact, max gap 1.108e+00 | IMPLEMENTED |
| A5 | Negative features break exactness | as A1 | R6 6/25 exact, max gap 9.195e−01 | IMPLEMENTED |
| A6 | Exactness is insensitive to n, k, depth, width, feature scale | `f5_stress_extremality.part_a_trial` | `runs/20260904-034701`: sum 400/400 live, 0 violations; max 286 live, 0 violations. Re-run `runs/20260905-122429` identical | IMPLEMENTED |
| A7 | max aggregation goes void as the mandatory subgraph densifies | `f5.part_b_live_frac` | same run: live fraction 1.00 → 0.00 as mandatory edges 0 → 23 | IMPLEMENTED |
| A8 | That failure mode does not reach RNA | `f5.part_c_rna` | same run: max mandatory degree 1, median density 0.28 % | IMPLEMENTED |
| A9 | Exactness survives residuals, depth 1–5, all monotone activations | `c6_closure.exact_frac` | `runs/20260904-034625`: 52/52 live exact, 0 violations, 2 voids. Re-run `runs/20260905-122446` identical | IMPLEMENTED |
| A10 | A non-monotone activation breaks exactness (discriminative control) | as A9 | same run: `sin` 0/15 exact in all 4 configs | IMPLEMENTED |
| A11 | Endpoint exactness ⟺ aggregator monotone under multiset inclusion, either direction | `aggregators.direction`, `c3.trial` | `runs/20260904-034519`: 7/7 agree; `min` maximum at the **empty** set | IMPLEMENTED *(empirical — see B1)* |
| A12 | k grows with length; 2^k intractable at 150 nt | `f2_lattice_size.main` → `ensemble.bpp_lattice` | `runs/20260904-020404`: median k 16.0 / 48.5 / 86.0 at 40 / 80 / 150 nt | IMPLEMENTED |
| A13 | Trained certifiable model reaches Spearman 0.592 held out | `train.fit`, `f6_trained.main` | `runs/20260904-013609`: ρ 0.5916, r 0.5564, n = 249 | IMPLEMENTED |
| A14 | Gene-disjoint Spearman 0.543 | `f6_trained` gene-disjoint block | same run: ρ 0.5432, r 0.5307, n = 290 | IMPLEMENTED |
| A15 | Published split shares all 30 genes | same | same run: `genes_shared: 30` | IMPLEMENTED |
| A16 | Sign constraint costs 0.002 Spearman (positional), 0.214 (generic) | `f8_baseline.evaluate` via `train.fit` | `runs/20260904-015628`, 5 seeds each | IMPLEMENTED |
| A17 | At floor 0 + canonical, bound holds on 100 % of sampled structures | `ensemble.sound_lattice`, `f7_soundness.main` | `runs/20260904-020253`: coverage 100.0 %, 0 / 20000 violations | IMPLEMENTED *(empirical — see B2)* |
| A18 | Soundness costs 29.25× at 60 nt | as A17 | same run, floor 0 canonical row | IMPLEMENTED |
| A19 | Floor 0.05 is not merely unsound in principle: 15 real structures exceeded it | as A17 | same run, floor 0.05 row: 15 / 20000 | IMPLEMENTED |
| A20 | 10^1183 structures bounded in 2 forward passes on real target sites | `f9_target_context.main`, `certify.certify_threshold` | `runs/20260904-020229`: 150 nt, k 3932, sound on every site. Re-run `runs/20260905-122512` identical | IMPLEMENTED |
| A21 | MFE understates the sampled ensemble max by 4.4–5.1 % | as A20 | same run, all three lengths | IMPLEMENTED |
| A22 | Maximal-element bound is near-exact (1.008–1.045) but costs M passes | `maximal.saturate`, `f10_maximal.main` | `runs/20260904-023117` | IMPLEMENTED |
| A23 | The maximal set does not converge at 150 nt (N^0.80) | `f10b_convergence.main` | `runs/20260904-023333`: exponents [0.754, 0.855, 0.902, 0.751] | IMPLEMENTED |
| A24 | Relaxation gap equals a walk-count ratio κ_L in the no-bias uniform-feature regime | `kappa.kappa_L`, `kappa.walk_count` | `runs/20260904-043204`: max relative error 2.121e−16, K1–K4 all pass | IMPLEMENTED |
| A25 | κ₂ bounds the measured slack and tracks its growth | `c8_kappa_gencode.main` | `runs/20260904-043511`: bounded at all 3 lengths; growth ×4.17, ×2.22 vs measured ×4.39, ×2.36 | IMPLEMENTED |
| A26 | The slack exponent is set by depth (0.912 per layer), not by RNA | `c5_slack_scaling.fit` | `runs/20260904-034559`: 0.912 at depths 1–4, R² 0.9996 | IMPLEMENTED |
| A27 | Slack is unbounded in problem size | `c5.slack` on matchings | same run: 9.0 (m=3) → 2290.3 (m=48) | IMPLEMENTED |
| A28 | Join-closed families have slack exactly 1.0 | `c5` part C, brute force | same run: `part_c_exact: true` | IMPLEMENTED |
| A29 | κ_L bounds the gap outside its derivation regime | `c9_kappa_boundary.run_case` | `runs/20260904-043246`: `kappa_bounds_all_cases: true` | IMPLEMENTED |
| A30 | Huesken table matches published statistics | `data.verify`, `data.load_dsir_split` | `make data` 2026-09-05: 2182/249, overlap 0, `verified: true` | IMPLEMENTED |
| A31 | Top-k rank certification is withdrawn | `certify.certify_topk` raises | `runs/20260904-015524` records `status: retired` | IMPLEMENTED |

---

## B. Implemented but weaker than the wording suggests

These are supported by evidence but the surrounding prose overstates the epistemic status.
Each needs a wording fix before submission, not a new experiment.

| # | Claim as written | What is actually established | Required correction |
|---|---|---|---|
| B1 | "Endpoint exactness holds **if and only if** the aggregator is monotone under multiset inclusion" | 7 aggregators × 25 seeded trials at n = 7, k = 8, depth 2, relu. **No proof.** The "only if" direction rests on three counterexample families | State as an empirically supported characterisation, or prove it. It is very likely provable in one direction by a two-line monotonicity argument; the converse needs care |
| B2 | The floor-0 certificate is "sound over the ensemble" | 0 violations in 20 000 *sampled* structures (A17), plus a monotonicity argument in the docstring | The containment argument (every valid structure's pairs ⊆ the canonical set) is sound and easy to prove. It is not written down as a proof anywhere. Do so |
| B3 | `exact_interval`, "THE THEOREM", "exact" throughout | float64 arithmetic, tolerance 1e-12, no interval or rational arithmetic, no directed rounding | Per the PI's standing instruction against unsupported "exact" claims: either add interval arithmetic, or qualify every use as *combinatorially exact, numerically float64* |
| B4 | "0 violations across 738 live trials" (report headline) | True, but every validation ran at k ≤ 11 while deployment runs at k = 408–3932 | State the validated range explicitly. Gap G-4 in the architecture audit |
| B5 | F5's independent reimplementation "shares no code with certmp" | True when run, but the file was never committed | Commit it, or drop the independence claim |

---

## C. Proposed — no code exists

All from `docs/stage1_novelty_audit.md`. Confirmed absent by the function inventory in
`docs/architecture_audit.md` §4.

| # | Proposed component | Nearest existing code | Status |
|---|---|---|---|
| C1 | Exact learned additive structured anchor (DP over the feasible family) | none — no DP in the repository | PROPOSED |
| C2 | Signed-weight nonlinear adjacency-dependent MPNN as the *certified* object | `MonoMPNN(nonneg=False)` exists but is **rejected** by the certificate gate | PROPOSED |
| C3 | Certified conditionable graph-family support oracle | `kappa.enumerate_family` (2^20 cap), `maximal.saturate` | PROPOSED |
| C4 | Edge-conditioned message-bound layer | none | PROPOSED |
| C5 | Learned allocation of a fixed support-query budget | none | PROPOSED |
| C6 | Certificate schema (assumptions, numerical status, hashes, bounds, witness, termination) | `provenance.new_run` records config + hashes but no certificate fields | PROPOSED |

---

## D. Unsupported — asserted but contradicted or superseded

| # | Location | Assertion | Why unsupported | Severity |
|---|---|---|---|---|
| D1 | `README.md:9` | Certificate bounds "predicted **off-target** risk" | The model is fit to Huesken **efficacy** — `%Inhibition` of the *intended* target. No off-target quantity is computed anywhere | **High** — misdescribes the application |
| D2 | `README.md:11` | Bound "loose by a measured median factor of 1.82" | 1.82× is the superseded **banded, unsound** lattice figure (F3b, untrained model). Sound figures are 29.25× and 222.8× | **High** |
| D3 | `README.md:18` | "threshold certificate + exact top-k rank stability" | `certify_topk` raises `NotImplementedError` | Medium |
| D4 | `README.md:2-4` | "provided the aggregation is sum or max" | Superseded by A11: logsumexp and min are also certifiable | Medium |
| D5 | `README.md:14-20` | module listing | Omits `aggregators.py`, `train.py`, `maximal.py`, `kappa.py` | Low |
| D6 | `RESULTS.md:398` | F10 100 nt covered closure "10^10 to 10^13" | Run record gives 10^10.235 – **10^12.478** | Medium — propagated into the shared report |
| D7 | `certmp/ensemble.py:5` | pairs above the band are "present in essentially every structure" | F3b measured median 96.6 %, worst 83.0 % of samples containing all mandatory pairs | Medium — this is the assumption that made the lower bound unsound |
| D8 | Report §"first"-adjacent phrasing | none currently uses the word "first" | — | none found; keep it that way |

**No claim in section A depends on any item in section D.** The unsupported items are all
documentation drift in `README.md` plus one rounding error, none of which feeds a computed
result.

---

## E. Ledger summary

| status | count |
|---|---|
| IMPLEMENTED and supported | 31 |
| IMPLEMENTED but overstated in prose | 5 |
| PROPOSED, no code | 6 |
| UNSUPPORTED | 7 |

The experimental core is in good order. The liabilities are concentrated in `README.md`,
which has not been updated since commit `45a72f0` and now contradicts the code in four
places, and in the word "exact", which is used in a combinatorial sense throughout while
the arithmetic is unverified float64.
