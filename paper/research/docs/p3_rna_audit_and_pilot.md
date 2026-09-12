> Publication copy of `docs/p3_rna_audit_and_pilot.md`. Mathematical and result content is retained; execution material is omitted. This dated record retains its original conclusions and limitations, including conclusions superseded by later stages. Original SHA-256: `0dc1e5b4ac4bf541d638d85fe35a4cf5d402f0e3cc46f37b1ecc81906f18134b`.

# P3 — RNA data audit and on-target efficacy pilot

Records: `runs/revision/rna_audit-20260906-003113.json`,
`runs/revision/rna_pilot-20260906-011411.json` (+ `.sha256`).
Configs: `configs/revision_rna_data.yaml`, `configs/revision_rna.yaml`.

**Scope: on-target efficacy and candidate ranking.** Not off-target safety, not
causal structure validation, and no thermodynamic claim for modified chemistry —
this dataset has none.

## The finding that changed the position

The previous report concluded that **six of eight arms were blocked** because the
reporter construct context was unrecoverable. That was premature. Before
declaring it lost, I mapped each guide's reverse complement into GENCODE isoforms
and took, per gene, the convex hull of the tiled site positions:

| Gene | Sites | Transcript nt | Tiled span | % of transcript |
|---|---:|---:|---:|---:|
| UBE2B | 79 | 1064 | 249–828 | 54.4% |
| UBE2E3 | 79 | 909 | 173–746 | 63.0% |
| UBE2D3 | 78 | 3729 | 303–643 | 9.1% |
| UBE2C | 76 | 641 | 249–489 | 37.4% |
| P2RX3 | 90 | 3792 | 237–1364 | 29.7% |

The siRNAs tile a **dense contiguous span** per gene, not the whole transcript —
the signature of an inserted fragment. That hull is a **lower bound on the
insert**, so a window lying wholly inside it has native flanking sequence that is
inside the insert:

| Window | Interior | Edge | Interior % |
|---:|---:|---:|---:|
| 50 nt | 1156 | 134 | **89.6%** |
| 100 nt | 1001 | 289 | 77.6% |
| 150 nt | 844 | 446 | 65.4% |

**Six of eight arms became admissible**, up from two — and admissibility is
decided by resolved inputs, not by dataset name.

**The assumption this rests on**, stated wherever it is used: the insert is a
contiguous native fragment. That is inferred from the tiling pattern, **not read
from a construct record**. The hull is a lower bound, so edge windows are
`unresolved_context` — not proven to be vector sequence. Vector backbone is not
recovered at all, and only 1290/2431 sites (53.1%) are locatable in GENCODE.

## What the audit preserved and refused

| Item | Status |
|---|---|
| Efficacy labels | **Unclipped.** 134 rows (5.51%) exceed 1.0, max 1.341, preserved |
| Assay identity | Primary source governs: **YFP reporter, H1299, 48 h**. The redistributed annotation's HeLa/luciferase columns are recorded and **contradicted**, never used as identity |
| F9 native windows | **Quarantined** — labelled structural examples only, supporting no measured reporter conclusion |
| Accession substitution | **Refused.** No plausible accession chosen, no construct invented |
| Isoform ambiguity | **906 of 1153** admissible rows match more than one isoform. Recorded per row, never folded into an exon union |

## Davis 2025 — still missing, and exactly what is needed

| | |
|---|---|
| Required file | `gkaf479_supplemental_files.zip` (contains **Supplementary Table S1**) |
| Route | Oxford Academic (*Nucleic Acids Research*) article page → Supplementary data |
| Place at | `data/davis2025/gkaf479_supplemental_files.zip` |
| Licence | CC-BY-4.0 — redistribution permitted once obtained |
| Why not automated | The PMC mirror serves a reCAPTCHA. **It was not bypassed and no challenge HTML was parsed.** |

Blocked until supplied: chemical scaffold covariates, chemical-transfer
evaluation, native-vs-reporter separation, dose/time covariates, replicate
structure, cross-dataset de-duplication. GEO **GSE231101** is 3P-seq
polyadenylation mapping and must never substitute for knockdown labels.

## Pilot result — the structural hypothesis is not supported

1,153 admissible interior windows at 50 nt, novel-sequence split (737 groups by
gene × 13-mer cluster, **group leak 0**), 5 seeds, metrics predeclared before any
arm ran.

| Arm | Spearman | sd | Δ vs sequence | P@20 | Enrichment |
|---|---:|---:|---:|---:|---:|
| matched_sequence_chemistry | 0.6595 | 0.0253 | — | 0.950 | 2.02 |
| target_sequence | 0.6595 | 0.0253 | +0.0000 | 0.950 | 2.02 |
| mfe_gnn | 0.6609 | 0.0296 | +0.0015 | 0.940 | 2.00 |
| sampled_ensemble_gnn | 0.6607 | 0.0294 | +0.0013 | 0.940 | 2.00 |
| anchor_A | 0.6609 | 0.0296 | +0.0015 | 0.940 | 2.00 |
| repaired_A_plus_B | 0.6607 | 0.0294 | +0.0013 | 0.940 | 2.00 |

Prevalence at the predeclared 0.70 activity threshold: **0.472**.

**Target structure adds +0.0015 Spearman, which is 5% of the seed noise**
(sd 0.028). Paired per-seed deltas are `[−0.0002, +0.0003, +0.0079, +0.0076,
−0.0083]` — sign changes across seeds. This is a null result, and now a *real*
one: the previous report could not run this test at all.

`target_sequence` equals `matched_sequence_chemistry` **exactly**, as it must —
the site is the reverse complement of the guide and carries identical
information. That is a consistency check, not a finding about target context.

## Prior sensitivity dwarfs the arm differences

| Prior temperature | Mean pair density |
|---|---:|
| 0.5 | 0.1728 |
| 1.0 | 0.1819 |
| 2.0 | 0.2001 |

Spread **0.0273** across priors — a *modelling* uncertainty, not sampling error,
and it is larger in relative terms than any difference between arms. Sampling
converges properly: standard error 0.168 at n=64 → 0.084 at n=256, halving as
√n predicts.

## Implementation limitation I must disclose

The arms named `mfe_gnn`, `anchor_A`, `sampled_ensemble_gnn` and
`repaired_A_plus_B` are implemented as **ridge regression on guide one-hot plus
four window-derived summary features** (pair density, mean span, mean marginal,
candidate density), with the sampled ensemble term added for the γ=1 arms. They
are **not** the full `CertifiedModel` graph network run over each window family.

The names therefore overclaim, and the correct reading is *structure-informed
feature* arms. Given P1 found a pairwise ridge beats the full residual outright,
using ridge here is defensible on cost grounds — but it means this pilot tests
**whether ensemble-derived window summaries carry signal**, not whether the full
architecture does. A full-model version is affordable and is the obvious next
step if the summaries had shown anything; they did not.

Label semantics are correct regardless: the assay outcome is fitted through
`E_p[f] = b + aᵀμ + γ·E_p[g]`, and no individual latent structure is ever given
the assay outcome.

## Where this leaves the RNA study

The honest position, with the structural claim separated from the sequence one:

- **Sequence-only prediction is solid**: Spearman 0.66, P@20 0.95, enrichment
  2.0 on a clean novel-sequence split with zero group leak.
- **Target structure adds nothing measurable** on the admissible data, at an
  effect 20× smaller than seed noise.
- **Chemical transfer and native-vs-reporter remain untestable** without Davis S1.
- Three caveats bound even the null: the contiguity assumption, 78% isoform
  ambiguity, and 53% site locatability.

This does not establish that target structure is irrelevant to siRNA efficacy.
It establishes that on 1,153 windows whose context is recoverable under a stated
assumption, ensemble-derived structural summaries do not improve on-target
efficacy prediction over guide sequence alone.
