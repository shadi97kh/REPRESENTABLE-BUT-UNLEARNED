from pathlib import Path
import shutil,json,hashlib,ast,subprocess
R=Path('/home/shadi/iclr2027');W=Path('/tmp/iclr2027-sirna-v3-publication');run=R/'runs/sirna_gnn_empirical/v3-20260914T013314Z';P=R/'papers/interaction_recoverability_iclr2027/v8';prefix='sirna_gnn_empirical/results/v3';selected={}
def copy(src,dest):
 dst=W/dest;dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(src,dst);selected[dest]=str(src.relative_to(R))
def write(dest,txt):
 dst=W/dest;dst.parent.mkdir(parents=True,exist_ok=True);dst.write_text(txt);selected[dest]='publication documentation'
names='activity analysis_supplement architecture audit_results b3_admission b3_evaluate b3_primary checks common engine evaluate figures legacy_analysis main_cell_provenance metrics pair pair_gradient_check parameter_audit protocol runner source_audit_finish source_checks'.split()
for n in names:
 src=R/f'sirna_gnn_empirical/v3/{n}.py';ast.parse(src.read_text());copy(src,f'sirna_gnn_empirical/v3/{n}.py')
for n in ['active_run.txt','run_all.sh']:copy(R/f'sirna_gnn_empirical/v3/{n}',f'sirna_gnn_empirical/v3/{n}')
# Retain every aggregate table; omit raw strand/endpoint/condition records.
excluded={'B2_chemistry_templates.csv','B2_duplex_states.csv','B2_measured_endpoint_mapping.csv','b3_primary_measurements.csv','b3_primary_states.csv'}
for src in (run/'tables').iterdir():
 if src.is_file() and src.name not in excluded:copy(src,f'{prefix}/tables/{src.name}')
for n in ['protocol.json','protocol.md','protocol_freeze.json','experiment_matrix.csv','activity_selection.json','pair_selection.json','environment.json','final_resource_accounting.json','resource_accounting.json','pair_gradient_checks.json','parameter_audit_summary.json','analysis_supplement_summary.json','main_cell_provenance.json','figure_provenance.json','table_provenance.json']:
 copy(run/n,f'{prefix}/{n}')
for n in ['S1_status.json','primary_release_discrepancies.json','reproduction_status.json','geometry_status.json','official_pin.json','megmod_unresolved_forward_symbols.json']:
 copy(run/'sources'/n,f'{prefix}/sources/{n}')
for n in ['admission_summary.json','evaluation_summary.json','protocol.json','source_holdout_checks.json']:
 copy(run/'b3_primary'/n,f'{prefix}/b3_primary/{n}')
for src in (P/'figures').iterdir():
 if src.suffix in ['.pdf','.png','.svg']:copy(src,f'sirna_gnn_empirical/presentation/v3/figures/{src.name}')
for src in (P.with_name('v7')/'figures').iterdir():
 if src.suffix in ['.png','.svg'] and not src.name.startswith('workflow.'):
  copy(src,f'sirna_gnn_empirical/presentation/v3/figures/{src.name}')
# Aggregate plotted marks only. Individual plotted measurements remain in the complete local archive.
for stem in ['main_support','main_activity','appendix_learning_curves','appendix_seed_risk','appendix_ranking_pools']:
 src=run/'figures'/f'{stem}_marks.csv'
 if src.exists():copy(src,f'sirna_gnn_empirical/presentation/v3/figure_data/{src.name}')
write('sirna_gnn_empirical/presentation/v3/README.md','''# Current workflow and result figures

These are the existing campaign v3 figures used by manuscript v8, **Do Graphs Learn the Chemistry? Controlled Evidence from siRNA Prediction**. No fitting or scientific reevaluation was performed to publish them.

![Architecture and separate supervision protocols](figures/workflow.png)

The supplied 3300 × 1950 PNG is preserved byte for byte. The standalone PDF embeds that image without a Figure caption; it is a raster PDF, not a vector drawing. The superseded SVG is not presented as its source. In the paper, Figure 1 is on page two with a separate caption below it.

The self transform is drawn with column vectors; the code uses the equivalent transposed row-vector convention. R0 uses eight typed transforms. R1 uses a shared directional backbone base and a residual gated by positive training occurrence. The no-message control substitutes receiver states while retaining degrees. Activity loss is equal-component standardized squared error. Pair supervision adds lambda times equal-component squared error of `(z1-z0) - (y1-y0)/s_train`, with shared weights; grouped development selects lambda. Outer held-out labels enter evaluation only; inner validation selects models. Exact equations are implemented in [architecture.py](../../v3/architecture.py) and [engine.py](../../v3/engine.py).

| Figure | Content |
|---|---|
| [Main support](figures/main_support.pdf) | Directed relation counts and training gradients |
| [Main activity](figures/main_activity.pdf) | Ten-seed errors, ensembles and calibration |
| [Main B2/ranking](figures/main_pair_ranking.pdf) | Tie-aware APP ranking and measured chemistry differences |
| `appendix_*.pdf` | Eight existing supporting figures: learning, calibration, seeds, ranking, B2, B3 and source reconciliation |

Each result figure has PDF, SVG and PNG exports. [Plotting code](../../v3/figures.py) is preserved as executed; it reads original run paths and copies figures into the historical local manuscript directory, which is not included in GitHub. It requires omitted row-level inputs for some panels and is not a standalone reproduction command for this aggregate checkout. Aggregate plotted marks are in `figure_data/`; full marks and predictions remain in the local scientific archive. Source paths/hashes are recorded in [figure provenance](../../results/v3/figure_provenance.json). Figure assets visualize recorded measurements/predictions; they are not new measurements.
''')
write('sirna_gnn_empirical/v3/README.md','''# Empirical campaign v3: original executed scientific code

The completed campaign ran 2,312 fits: 1,776 development and 536 final, with ten final seeds for stochastic methods. It recorded 1,407,967 executed neural optimizer updates. No core fits failed; unfavorable finite seeds are retained. Exact ENsiRNA-mod/MEG-mod supervised reproductions remain blocked, with zero such fits. The separate primary B3 evaluation reused 62 source-held-out checkpoints and performed no fitting.

The source files here are byte-preserved copies of the executed code. `run_all.sh` uses the original `/opt/miniforge3/bin/python`; `common.py` resolves the original run via `active_run.txt` and reads the historical v1/v2 inputs. Frozen [protocol/configurations](../results/v3/protocol.json), [selection scores](../results/v3/tables/selection_scores.csv), [fit records](../results/v3/tables/all_fit_summary.csv) and [software versions](../results/v3/environment.json) are included.

This focused Git checkout contains code, aggregate results and presentation assets, not the measured rows, graph tensors, all predictions or checkpoint binaries required to rerun the campaign. The original local run remains `runs/sirna_gnn_empirical/v3-20260914T013314Z/`, with the complete `paper_artifacts/code_results.zip`. Missing local files/completion receipts in GitHub do not authorize or trigger a rerun.

The original resumable entry point was `bash sirna_gnn_empirical/v3/run_all.sh`. Its scientific stage order is `legacy_analysis source_checks protocol checks activity pair b3_evaluate evaluate analysis_supplement parameter_audit figures`; default `paper`/`deliver` stages depend on local-only manuscript files and are excluded from this publication. Acquired primary-source follow-up is specified by `b3_primary.py`, `b3_admission.py` and the retained source checks. Original completion receipts and command logs stay in the complete local campaign. No scientific module was imported or executed while preparing this Git contribution. For an intentional reproduction, first provide the documented inputs, use a fresh run directory and preserve the frozen split/selection protocol; do not overwrite the completed run.

The GNN is conventional, with no new layer, no five-law guarantee and no established general advantage over the strongest fitted baseline. APP/S7 are retrospective follow-up. B2 has not established sequence-specific measured chemistry prediction. See [current results and limitations](../README.md).
''')
old=subprocess.check_output(['git','show','origin/main:sirna_gnn_empirical/README.md'],cwd=W,text=True)
current='''# Measured siRNA GNN: code and controlled comparisons

The latest completed study is **empirical campaign v3**, presented as *Do Graphs Learn the Chemistry? Controlled Evidence from siRNA Prediction*. The conventional chemistry-aware GNN was actually trained against published measured outcomes. The campaign tests relation-support correction crossed with message passing, ten final seeds, grouped selection, matched controls, ranking and measured chemistry differences.

| v3 row-weighted MSE | Grouped ENsiRNA | APP | Davis S7 |
|---|---:|---:|---:|
| Corrected GNN | .071175 | .162475 | .109303 |
| Corrected no-message | .071237 | .163043 | .104709 |
| Original no-message | .071237 | .162914 | .105706 |
| Chemistry tree | .065648 | .189569 | .091901 |
| Token CNN | .071749 | .164365 | .101685 |
| Training row mean | .072372 | .175983 | .148803 |

No general GNN improvement is established. On APP, GNN minus tree MSE is −.027094 [−.041964, −.010482], but GNN minus the strongest observed fitted baseline (original no-message) is −.000439 [−.001856, .001151]. GNN APP R-squared is −.001360 and prediction SD .007201. These conditional component intervals hold fitted predictions fixed; they do not include training-seed variability or certify new-study generalization. Grouped ENsiRNA remains a released-label benchmark; APP and S7 are retrospective follow-up.

B2 has 156 pairs in 12 sequence components. Pair-GNN MSE .222115 compares with matched no-message .222013, pair ridge .219297, training-fold mean .263770 and zero .280804. Improvement over the fold mean remains unresolved; sign accuracy matches the majority rule. Sequence-specific measured chemistry prediction was not established. The joint four-position bundle cannot isolate GNA.

Primary Bramsen B3 evidence contains 165 measured states and 140 four-condition rectangles on one sequence background. Fixed-checkpoint GNN interaction MSE .068335 does not improve on zero .068205. APP still has no complete rectangles. Common-reference dependence, missing endpoint replicate/covariance information and primary/released-label discrepancies remain limitations: 299 of 1,604 compared values differ by more than .05, with maximum 1.0852. Acquired S1 remains quarantined for unresolved orientation.

## Current artifacts

- [Executed scientific code and reproduction scope](v3/README.md).
- [Complete aggregate comparison/fit tables](results/v3/tables/) and [frozen protocol](results/v3/protocol.json).
- [Current workflow and all eleven result figures](presentation/v3/README.md).
- [Published-method blockers](results/v3/sources/reproduction_status.json), [source reconciliation](results/v3/sources/primary_release_discrepancies.json), and [resource accounting](results/v3/final_resource_accounting.json).
- [Current publication manifest](../docs/sirna_gnn_empirical/publication/v3/files.json).

The core campaign executed 2,312 fits (1,776 development, 536 final), with ten final seeds for stochastic methods and 1,407,967 neural optimizer updates. Exact published-method fits remain zero. B3 evaluation reused source-held-out checkpoints without fitting. No new wet-lab validation, causal mechanism, architectural novelty or five-law theorem coverage is claimed. Backend no-go and the historical reference-efficiency question are unchanged.

Run `python sirna_gnn_empirical/verify_published.py` for read-only verification of the current published files. This publication includes code, aggregate results, workflow and result figures. Full measured rows, per-observation predictions, checkpoint binaries, original command logs and complete archives stay preserved locally. Paper TeX/manuscript PDFs and private authoring material are not pushed. Reproduction limitations are explicit in the code README.

## Preserved historical v1/v2 publication

The material below describes those earlier campaigns only; its B3 status predates the primary-panel follow-up above. Historical metrics and negative results are unchanged.

'''
write('sirna_gnn_empirical/README.md',current+old.replace('# Measured siRNA GNN: code and aggregate results','### Original publication',1))
root=subprocess.check_output(['git','show','origin/main:README.md'],cwd=W,text=True);intro='''## Current empirical siRNA study

**Do Graphs Learn the Chemistry? Controlled Evidence from siRNA Prediction** evaluates a conventional chemistry-aware GNN against published measured outcomes. The latest campaign completed 2,312 fits with ten final seeds for stochastic methods. It does not establish a general GNN advantage over the strongest fitted baseline or sequence-specific measured chemistry prediction; unfavorable results and source limitations remain explicit. APP/S7 evaluation is retrospective, and no new wet-lab validation was performed.

[Current results, code and limitations](sirna_gnn_empirical/README.md) · [Comparison tables](sirna_gnn_empirical/results/v3/tables/) · [Workflow and result figures](sirna_gnn_empirical/presentation/v3/README.md)

![Chemistry-aware siRNA architecture and separate training protocols](sirna_gnn_empirical/presentation/v3/figures/workflow.png)

The supplied workflow is the paper's page-two figure. Its notation, controlled support/no-message variants and pair-loss normalization are documented with the [figure assets](sirna_gnn_empirical/presentation/v3/README.md). The empirical GNN is separate from the theoretical estimator described below.

'''
root=root.replace('## Current research result',intro+'## Separate theoretical research result',1)
write('README.md',root)
write('sirna_gnn_empirical/presentation/README.md','''# siRNA presentation assets

The latest workflow and all eleven result figures are in [v3/](v3/README.md), matching empirical campaign v3 and manuscript v8. The supplied workflow includes separate activity and pair-supervised training lanes.

The earlier v2 assets remain preserved in this directory. Their workflow illustrates the original activity-only strip; it is superseded for the current paper. Earlier `render_main_figures.py` and `figure_data/` describe the existing v2 results, not campaign v3. No paper TeX or manuscript PDF is published here.
''')
write('sirna_gnn_empirical/verify_published.py','''"""Read-only SHA-256 verification of the current scientific publication."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent.parent
manifest={}
# Later manifests override only explicitly revised documentation/assets.
for version in ['v1','v3']:
    manifest.update(json.loads((ROOT/f'docs/sirna_gnn_empirical/publication/{version}/files.json').read_text()))
for name,expected in manifest.items():
    p=ROOT/name
    if not p.is_file():raise SystemExit(f'Missing published artifact: {name}')
    h=hashlib.sha256()
    with p.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024),b''):h.update(block)
    if p.stat().st_size!=expected['bytes'] or h.hexdigest()!=expected['sha256']:
        raise SystemExit(f'Published artifact integrity mismatch: {name}')
print(f'Verified {len(manifest)} current published scientific artifacts; no training or inference.')
''')
write('docs/sirna_gnn_empirical/publication/v3/README.md','''# Campaign v3 scientific publication

Explicit user authorization: push necessary code to GitHub, update the README and include the supplied workflow and result figures. Destination: `shadi97kh/iclr2027`, existing `main` branch. Prepared from the remote main history in an isolated checkout, preserving the original local branch, all campaigns and manuscript versions. No force push is intended.

`files.json` hashes the added/revised files; the verifier overlays these entries on the historical v1 manifest. `source_mapping.json` identifies byte-preserved local sources. Published material comprises scientific scripts, configurations, aggregate comparison and per-fit summaries, resource/source-status records, the supplied caption-free workflow PNG/PDF and existing result figures. Figure assets include the requested scientific plots, while underlying per-observation prediction files, raw datasets and checkpoint binaries remain local. Paper TeX/PDFs, review archives, task instructions and dependency caches are excluded.

No fitting, inference, new scientific analysis or historical ledger write occurred during publication. Compilation, copying, hashing, Git operations and network transfer incur nonzero administrative cost, not a separately metered scientific resource charge.
''')
write('docs/sirna_gnn_empirical/publication/v3/source_mapping.json',json.dumps(selected,indent=2)+'\n')
manifest={name:dict(bytes=(W/name).stat().st_size,sha256=hashlib.sha256((W/name).read_bytes()).hexdigest()) for name in selected}
dest='docs/sirna_gnn_empirical/publication/v3/files.json';(W/dest).write_text(json.dumps(manifest,indent=2)+'\n')
(R/'docs/sirna_gnn_empirical/publication/v3/allowlist.json').write_text(json.dumps(list(selected)+[dest],indent=2)+'\n')
print('selected',len(selected)+1,'bytes',sum(x['bytes'] for x in manifest.values()))
print('source code byte-preserved',len(names))
