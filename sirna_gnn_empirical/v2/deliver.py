"""Administrative authoring and packaging of already-completed results."""
from core import *
import re,shutil,zipfile,subprocess
D=ROOT/'docs/sirna_gnn_empirical/v2';P=ROOT/'papers/interaction_recoverability_iclr2027/v3';W=P/'workflow_handoff';W.mkdir(exist_ok=True)
spec=dict(title='Executed targeted siRNA support and supervision follow-up',slot=dict(width_inches=5.5,height_inches=1.7,caption_below=True),status='Final illustration not drawn',lanes=[dict(name='Frozen v1 evidence',nodes=['4766 measured observations','Original grouped split and trained checkpoints','Support/gradient audit','Four frozen inference diagnostics'],color_role='preserved historical'),dict(name='Targeted B1 refits',nodes=['Supported representation x supported context (2x2)','Same original split; two-config selection; three seeds','Matched tree/token/no-message/exact chemistry controls','Retrospective APP prediction and ranking'],color_role='new executed'),dict(name='Grouped ENsiRNA',nodes=['2927 rows; ten study-linked components','Five outer folds; grouped inner selection','All four corrected models; all outer groups reported','Inner-selected all-ENsiRNA deployment fits'],color_role='new executed'),dict(name='Supervised B2',nodes=['No outside-APP pairs meet full comparability','156 APP pairs / 26 backgrounds / 12 sequence components','Four outer folds; paired endpoints remain together','Activity-only versus pair loss and fold-trained controls','No clear sequence-specific advantage over fold mean'],color_role='exploratory within cohort'),dict(name='Separate source',nodes=['Davis S1 remains quarantined','118 literal S7 duplexes pass structural and overlap checks','Frozen ENsiRNA deployment checkpoints','Qualified native-outcome prediction; dose/time masked'],color_role='qualified separate source')],blocked_nodes=['B3: zero complete measured reference/A/B/A+B rectangles','Five-law applicability: no verified private anchors or latent law','No new wet-lab, causal-mechanism or clinical validation'],required_distinctions=['APP follow-up is retrospective, not independent blind external evaluation','S7 is separately sourced observational evidence, not prospective validation','Fixed trained GNN contrasts are not measured biological effects','The empirical-quantile estimator and neural-profile/Riesz pilot are different procedures'],edges=['source checks -> label-independent grouping -> fitting/selection -> prediction commitment -> measured evaluation','grouped inner development -> frozen deployment -> S7 inference','B2 pair groups -> grouped training and validation -> outer-pair evaluation'],do_not_draw=['invented biological validation','B3 performance curve','private-anchor calibration that was not observed','new GNN novelty badge'],source_paths=['protocol.json','support_checks.json','factorial_fits.json','grouped_fits.json','b2_admission.json','b2_fits.json','external_admission.json','external_freeze.json','prediction_commit.json','external_prediction_commit.json'])
write(W/'workflow_spec.json',spec)
(W/'workflow_spec.md').write_text('''# Scientific workflow handoff

Create one legible vector workflow for a 5.5-inch-wide, 1.7-inch-high reserved slot. Do not put the caption inside the artwork. The current manuscript reserves the slot and does not contain a final illustration.

Use the five executed lanes in workflow_spec.json. Distinguish preserved v1 checkpoints from new B1 refits, grouped ENsiRNA evaluation, exploratory APP B2 supervision and qualified Davis S7 prediction. Shared arrows may connect measured-source verification and grouping, but never connect held-out outcomes back into selection. Show a visibly blocked B3 branch because no four-condition labels exist. No anchor calibration or wet-lab validation occurred.

Keep labels concise, preserve the exact cohort/group counts, and use consistent model/comparison terminology. The controlled architecture is conventional. Include no drawing that suggests mechanistic causality, clinical efficacy, new architecture or application of the five-law theorem to this GNN. APP follow-up and S7 must remain distinct. The JSON supplies the full source-grounded content and allowed relationships.
''')
(W/'workflow_source_excerpts.md').write_text('\n\n'.join('## '+name+'\n\n```json\n'+(RUN/name).read_text()+'\n```' for name in ['support_checks.json','b2_admission.json','external_admission.json','deployment_selection.json']))
(W/'workflow_caption.tex').write_text(r'\caption{Executed scientific workflow. Preserved measured sources and checkpoints support a training-support audit; new controlled fits use grouped development, and the APP bundle task uses separate within-cohort supervision. Davis S7 uses frozen ENsiRNA checkpoints. APP follow-up is retrospective, S7 is qualified separately sourced observational evidence, and no complete B3 or new wet-lab experiment is available.}\label{fig:workflow}'+'\n')
# Primary citation audits are inherited explicitly, not represented as 68 newly opened papers.
shutil.copy2(ROOT/'docs/sirna_gnn_empirical/v1/citation_audit.md',D/'inherited_empirical_citation_audit.md');shutil.copy2(ROOT/'papers/interaction_recoverability_iclr2027/v1/private/citation_audit.md',D/'inherited_theory_citation_audit.md')
(D/'citation_audit.md').write_text('''The compiled bibliography contains 68 entries, including the preserved corrections and theoretical references. Exact inherited primary-source checks are copied alongside this audit. This follow-up did not repeat full-text verification of every historical paper; it retains their checked scope and does not import new GNN guarantees.

Davis 2025: publisher full article and the original supplemental workbook were inspected again. S1, S6 and S7 were checked separately, preserving literal strands and all source rows. S1 remains quarantined; all 118 S7 rows pass the specified source-derived duplex checks. Dose/time are not invented for S7. Primary source: https://academic.oup.com/nar/article/53/12/gkaf479/8171869 .

Shmushkovich 2018: publisher HTML, Europe PMC primary XML and publisher supplementary PDF were acquired and hashed. They identify a separate 356-compound/17-gene study with a compatible endpoint design; missing per-compound tabular/modification details prevent admission. Primary source: https://academic.oup.com/nar/article/46/20/10905/5085976 . This supplementary source is linked in text; it is not a fitted benchmark.

ENsiRNA and MEG-mod remain prior chemistry-aware graph precedents. The historical MEG-mod access is abstract-level; no precise layer, loss or numerical-performance comparison is newly asserted. No uncertain pretrained checkpoint is scored as a clean matched baseline. All 68 compiled citation keys resolve; the official bibliography style and conference style are preserved.
''')
(D/'claim_audit.md').write_text('''| Claim | Evidence | Qualification |
|---|---|---|
| Specific original weights had no data support | relation_support.csv, feature_support.csv, context_support.csv, data_gradients.csv, frozen_replay.json | Full training-data gradients, not aggregate-only diagnostics; AdamW shrinkage is distinct. |
| Frozen fallbacks do not repair APP errors | B1_comparison.csv phase=frozen, all arms/seeds | Modified inference on old checkpoints; original outcomes remain unchanged. |
| Corrected GNN does not beat corrected tree on APP | B1_comparison.csv, paired_comparisons.csv | Retrospective original-split comparison, conditional sequence-component interval. |
| Grouped GNN has similar pooled error but worse equal-component error | grouped_comparison.csv and every component result | Ten recorded components, not ten independent laboratories; old internal outcomes known. |
| Pair supervision learns a common negative bundle direction | B2_supervised_comparison.csv and all fold predictions | Same direction correctness as fold majority; no clear gain over fold mean; no isolated GNA effect. |
| Davis S7 is usable for qualified native-outcome prediction | source_audit/davis_orientation_rows.csv, external_overlap_checks.csv, external_freeze.json | 118 literal duplexes; 76 sequence components; masked unresolved dose/time; public-source observational evaluation. S1 remains quarantined. |
| B3 remains incomplete | eligible_rectangles.jsonl in preserved v1; b2_admission.json | No measured four-condition rectangle; no synthetic labels. |
| Conventional corrected model, no architectural novelty | core.py and complete fit checkpoints | The empirical-quantile estimator and neural-profile/Riesz pilot remain separate. |
| Calibrated theorem scope retained | Full theoretical manuscript/proofs copied with namespaced references | Pointwise attainment is not efficient attainment or finite-risk upper bound; fixed hard profiles could be disclosed. Backend no-go and historical reference-efficiency question remain. |

The original seven-model negative table remains in the main text. The accepted review amendments—support failure, unsupported context/chemistry, correct B2 coverage, selected/total updates, qualified constants, coherent B2 figures and captions below artifacts—are integrated into fresh manuscript v3. No reviewed source archive was supplied; the user-provided review PDF and audit were preserved, and their substantive amendments were incorporated into the new source.
''')
(D/'appendix_audit.md').write_text('''The unrestricted appendix contains the complete current protocol and result tables, then the entire preserved empirical manuscript/evidence, then the entire calibrated theoretical manuscript and its dependency-complete appendices. Every main proposition, scientific equation, figure and table has a roadmap entry. All historical internal references were namespaced and all dependent input files copied. The full theory includes its assumptions, normalization, commutator bounds, global calibration, empirical-quantile rules, influence components, common-reference covariance, simultaneous remainder schedule, DQM/local-path arguments, saturation bounds, nonlinear testing path, product affinities and full-range polynomial-count corollary. No mathematical derivation terminates at a repository link.

All current model/input/optimization/loss/grouping/metric/uncertainty conventions are written out. Every model/arm/seed, every engine fit, all outer components and all APP ensemble assay pools are printed. Raw row identities, original strings, full checkpoint histories and prediction arrays accompany the scientific archive. Historical content is explicitly marked historical; current source-specific S7 admission does not retroactively rewrite S1 quarantine or v1 conclusions.

Captions are below all figures and tables; longtable captions follow the final segment. The official fonts, margins and style file were not reduced. Main pages are exactly 1–9; page boundaries and total appendix length are recorded in build_audit.json. Workflow artwork remains deliberately reserved, with four source-grounded handoff files.
''')
cmd=[json.loads(x) for x in (RUN/'commands.jsonl').read_text().splitlines()];fits=[json.loads(p.read_text()) for p in (RUN/'fits').rglob('fit.json')]
resource=dict(engine_fits=len(fits),direct_pair_ridge_fits=4,checkpoint_recorded_optimizer_updates=sum(f['optimizer_updates'] for f in fits),profile_optimizer_updates=20,possible_uncheckpointed_interrupted_updates=[0,21],selected_state_updates=sum(f['selected_updates'] for f in fits),recorded_child_cpu_s=sum(x['child_cpu_s'] for x in cmd),recorded_aggregate_stage_wall_s=sum(x['wall_s'] for x in cmd),maximum_child_rss_mib=max(x['max_child_rss_kib'] for x in cmd)/1024,administrative_cost='Additional nonzero unmetered direct inspection, source acquisition, authoring, preflight builds, rendering, hashing and packaging. Runner paper/source-check/integrity costs are included above; direct administrative work is not.',original_files_verified=5402,no_historical_allowance_or_ledger_used=True)
write(D/'resource_accounting.json',resource)
(D/'executed_commands.md').write_text('# Exact executed stage commands\n\nCommands below are verbatim runner-expanded commands; full logs and failed attempts remain in the run. No prepared-only fit is presented as executed.\n\n'+''.join('```sh\n'+' '.join(x['command'])+'\n```\n\nExit '+str(x['exit_code'])+'; wall '+str(round(x['wall_s'],3))+' s; child CPU '+str(round(x['child_cpu_s'],3))+' s.\n\n' for x in cmd)+'''Additional direct administrative commands included `source_checks.py`, `prepare_paper.py`, `paper_tables.py`, `tectonic --keep-logs --keep-intermediates .../v3/main.tex`, `pdftotext`, `pdftoppm`, source HTTP acquisition and preservation hashing. Source checks were subsequently recorded as their own runner stage. Their earlier direct administrative cost is not falsely included in the stage ledger.

Resumable entry point: `bash sirna_gnn_empirical/v2/run_all.sh`. In a scientific-only archive, pass `audit source_checks plan factorial grouped b2 external evaluate figures`; paper/integrity stages additionally require the manuscript/full original workspace. New-machine fitting uses the same scripts with a fresh run path and the preserved data; the included completed manifests verify current outputs without training again.
''')
# All primary tables as readable Markdown companions.
def md(df):
 df=df.copy();return df.to_markdown(index=False,floatfmt='.6f')
for file in ['B1_comparison.csv','grouped_comparison.csv','B2_supervised_comparison.csv','external_comparison.csv','paired_comparisons.csv']:
 df=pd.read_csv(RUN/file)
 if 'seed' in df:df=df[df.seed=='ensemble']
 (D/(Path(file).stem+'.md')).write_text(md(df)+'\n')
summary='''The targeted follow-up ran in a fresh v2 campaign. All 5,402 preexisting files—including the completed v1 campaign, original manuscripts, proofs, frozen pilot, ledgers and review packets—retain their hashes.

The support audit confirms 8,192 original PO/PS message parameters with zero full-training-data gradient. Frozen inference fallbacks worsen APP MSE slightly. Corrected representation/context refits do not establish a GNN advantage: GNN APP MSE **0.185125**, corrected tree **0.162656**, difference **+0.022469 [0.004601, 0.038227]** conditional on this cohort. Original GNN/tree results **0.174726/0.168597** remain unchanged.

Five grouped ENsiRNA folds evaluate all 2,927 observations/ten components. GNN/tree pooled MSEs are **0.067062/0.067574**; equal-component MSEs are **0.063994/0.060579**. The equal-component difference interval **[-0.009647, 0.016175]** establishes no consistent GNN improvement. Every source group and seed is reported.

Exploratory within-APP supervised B2 used 156 pairs, 26 exact backgrounds and 12 sequence components. Pair GNN versus fold-mean pair MSE is **0.239508/0.237508**; equal-component MSE is **0.318228/0.326584**, difference **-0.008356 [-0.046418, 0.019650]**. All 137 measured non-ties receive negative calls, with **111 correct**, matching the training-fold majority. Learning the common harmful bundle does not establish sequence-specific chemistry-effect prediction. B2 jointly changes four positions and cannot isolate GNA. Replicate/control covariance remains unresolved; no SE is fabricated.

Davis S1 remains quarantined. All 118 literal S7 duplexes pass source-derived complementarity checks and have no detected 13-mer overlap with the prior data; they form 76 sequence components. Frozen all-ENsiRNA models were evaluated on these additional native-outcome measurements: GNN **0.091278**, tree **0.099535**, token CNN **0.089326**, no-message **0.109195** MSE. Token CNN is descriptively best. S7 is separately sourced observational evidence, with unresolved dose/time masked; it is not prospective blind or new wet-lab validation. Another separate primary study (Shmushkovich 2018) was identified, but its acquired supplementary PDF does not expose the referenced full per-compound chemistry table.

The learned architecture remains conventional. No measured B3 rectangle, new biological mechanism, clinical benefit or GNN theorem guarantee is established. The empirical-quantile estimator, earlier neural-profile/Riesz pilot and current GNN remain distinct. Backend no-go and the historical reference-efficiency question remain unchanged.

The full manuscript has **main 1–9**, excluded statements **10**, references **11–15**, appendix **16–182 (167 pages)**. It contains four main figure slots, nine main tables and 68 bibliography entries. All captions are below their artifacts, including multipage tables. Complete current methods/results and preserved historical empirical and theoretical material are included. The workflow is specified, not drawn.

Executed: **210 engine fits plus four direct pair-ridge fits**; **39,306 checkpoint-recorded optimizer updates plus 20 profile updates**, with **0–21 possible additional uncheckpointed updates** from the interrupted worker. Selected-state update counts sum to 14,654 across distinct fits. Resource totals and the failed attempt are in resource_accounting.json; direct administrative work has additional unmetered nonzero cost. No historical resource wrapper or ledger was used.

Artifacts:
- [Full PDF](../../../papers/interaction_recoverability_iclr2027/v3/main.pdf)
- [Results-focused PDF](../../../papers/interaction_recoverability_iclr2027/v3/results_extract.pdf)
- [Editable manuscript source](../../../papers/interaction_recoverability_iclr2027/v3/manuscript_source.zip)
- [Workflow handoff](../../../papers/interaction_recoverability_iclr2027/v3/workflow_handoff.zip)
- [Scientific code/results/checkpoints](../../../runs/sirna_gnn_empirical/v2-20260913T220847Z/code_results.zip)
- [Exact commands](executed_commands.md), [claim audit](claim_audit.md), [citation audit](citation_audit.md), [appendix audit](appendix_audit.md), [build audit](build_audit.json), [integrity audit](integrity_audit.json).

No repository push, publication or external contact was performed. The latest follow-up explicitly prohibits publishing; the earlier unresolved push restriction was not retried or bypassed.
'''
build=json.loads((D/'build_audit.json').read_text());summary=summary.replace('16–182 (167 pages)',f"{build['appendix'][0]}–{build['appendix'][1]} ({build['appendix_pages']} pages)");(D/'execution_summary.md').write_text(summary)
# Keep a self-contained code snapshot and explicit source identities in the scientific run.
source=RUN/'completed_source';source.mkdir(exist_ok=True)
for p in Path(__file__).parent.glob('*'):
 if p.is_file():shutil.copy2(p,source/p.name)
write(RUN/'completed_source_manifest.json',{p.name:sha(p) for p in source.iterdir() if p.is_file()})
# Editable manuscript archive contains every TeX dependency; no machine-local private logs.
with zipfile.ZipFile(P/'manuscript_source.zip','w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
 for p in sorted(P.rglob('*')):
  if not p.is_file() or 'private' in p.relative_to(P).parts:continue
  if p.suffix in ['.tex','.bib','.bst','.sty','.md','.json','.sh'] or ('figures' in p.relative_to(P).parts and p.suffix in ['.pdf','.png','.jpg']):z.write(p,p.relative_to(P))
with zipfile.ZipFile(P/'workflow_handoff.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(W.iterdir()):z.write(p,p.name)
with zipfile.ZipFile(P/'completed_audits.zip','w',zipfile.ZIP_DEFLATED) as z:
 for p in sorted(D.iterdir()):
  if p.is_file():z.write(p,p.name)
print('Authored documentation and manuscript/workflow archives')
