"""Aggregate already frozen pilot artifacts; never fit or simulate data."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path
import numpy as np


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root',required=True)
    parser.add_argument('phase',choices=['report'])
    args=parser.parse_args();run=Path(args.run_root)
    docs=Path(__file__).resolve().parent
    summary=json.loads((run/'evaluation/summary.json').read_text())
    rows=list(csv.DictReader((run/'evaluation/raw_results.csv').open()))
    aggregates=summary['aggregate'];methods=list(aggregates)
    costs={m:float(np.mean([float(r['method_cpu_s_fully_charged']) for r in rows if r['method']==m])) for m in methods}
    errors={m:{r['id']:float(r['squared_error']) for r in rows if r['method']==m and r['squared_error']} for m in methods}
    paired={m:{'adaptive_lower_error_datasets':sum(errors['neural_adaptive'][i]<v for i,v in errors[m].items()),
               'comparison_datasets':len(errors[m]),'relative_mse_reduction':1-aggregates['neural_adaptive']['mse']/aggregates[m]['mse']} for m in methods if m!='neural_adaptive'}
    geos=[];solves=[];fits=[];paths=[];trajectories=[];lambdas=Counter();dimensions=Counter();preds=[]
    for folder in sorted((run/'predictions').iterdir()):
        preds.append(json.loads((folder/'predictions.json').read_text()))
        for file in sorted(folder.glob('fold-*.json')):
            record=json.loads(file.read_text())
            for kind,model in record['models'].items():
                geos.append(model['geometry_error']);fits.append(model['fit'])
                solves.extend(model['methods'].values())
                if kind=='neural':
                    adaptive=model['methods']['adaptive'];lambdas[str(adaptive['lambda_'])]+=1;dimensions[str(adaptive['dimension'])]+=1
                    paths.extend(adaptive['finite_path']['rows'])
                    trajectories.extend(x for x in adaptive['trajectory'] if x['accepted'])
    diagnostics={
        'rows':len(rows),'datasets':len(preds),'rotation_models':len(fits),'selected_corrections':len(solves),
        'mean_method_cpu_s_fully_charged':costs,'paired_comparisons':paired,
        'fitter_total_cpu_s':sum(p['total_cpu_s'] for p in preds),
        'fitter_total_wall_s':sum(p['total_wall_s'] for p in preds),
        'max_score_gram_refinement_difference':max(g['gram_difference_norm'] for g in geos),
        'max_target_gradient_refinement_difference':max(g['gradient_difference_norm'] for g in geos),
        'max_fitted_score_center_abs':max(g['centering_max_abs'] for g in geos),
        'max_source_crps_refinement_difference':max(f['source_quadrature_difference'] for f in fits),
        'max_selected_solver_relative_residual':max(s['solve_relative_residual'] for s in solves),
        'max_selected_regularized_condition':max(s['regularized_condition'] for s in solves),
        'adaptive_lambda_counts':dict(lambdas),'adaptive_dimension_counts':dict(dimensions),
        'finite_paths':len(paths),'source_compatible_finite_paths':sum(p['source_compatible'] for p in paths),
        'max_abs_finite_path_B':max(abs(p['B']) for p in paths),
        'max_abs_nonlinear_remainder':max(abs(p['nonlinear_remainder']) for p in paths),
        'accepted_enrichment_directions':len(trajectories),
        'witness_after_to_before_median':float(np.median([t['same_direction_residual_after']/t['witness'] for t in trajectories])),
        'max_truth_target_refinement_difference':max(json.loads(p.read_text())['target_quadrature_difference'] for p in (run/'private').glob('d*.json')),
        'numerical_status':'Refinement differences and finite witnesses, not rigorous enclosures or whole-class upper certificates',
        'figures_visually_inspected':['figure1_interaction_error.png','figure2_ablation_cost.png','figure3_witness_nonlinearity.png'],
        'image_viewing_note':'Default view_image hit filesystem sandbox launcher failure; original PNGs were read via budgeted base64 and displayed without editing.'}
    with (run/'evaluation/diagnostics.json').open('x') as f:json.dump(diagnostics,f,indent=2,allow_nan=False)
    lines=['# Actual pilot results','', '**Decision: incremental-known.** All 24 frozen datasets and all ten methods completed, yielding 240 paired prediction rows. The adaptive method reduced aggregate MSE by 3.5569% against the strongest complete control, below the preregistered 10% screening threshold. Three replicates per cell support a descriptive comparison only.','',
           'The proposed learned nonlinear theorem remains open. A small observed gain does not establish a new statistical method, calibrated uncertainty or a general computational advantage. No additional fits were run after final scoring.','',
           '## All methods','', '| Method | Aggregate MSE | Complete datasets | Mean CPU seconds, full shared cost |','|---|---:|---:|---:|']
    for m in sorted(methods,key=lambda m:aggregates[m]['mse']):
        lines.append(f"| {m} | {aggregates[m]['mse']:.9f} | {aggregates[m]['successful_datasets']}/24 | {costs[m]:.6f} |")
    strongest=summary['strongest_complete_control'];pair=paired[strongest]
    lines+=['',f"Adaptive had lower squared error than `{strongest}` on {pair['adaptive_lower_error_datasets']}/24 individual datasets. Against the same-fit fixed correction its aggregate reduction was {100*paired['neural_fixed']['relative_mse_reduction']:.4f}%. These counts are paired descriptive outcomes, without a significance or rate claim.", '',
            'The direct-target set control found feasible endpoints on all 24 datasets. Its endpoint search is still approximate and restricted to the basis sieve: it is neither an optimal original-class comparator nor a certified confidence interval. No control received true profiles, coefficients, latent observations or joint outcomes.', '',
            '## Cell-level comparison','', '| Family | delta | n per stratum | Adaptive MSE | Rich neural control MSE | Fixed neural MSE | Replicates |','|---|---:|---:|---:|---:|---:|---:|']
    cells=summary['per_cell']
    for family in ['oscillatory','localized']:
        for delta in [.1,.03]:
            for n in [256,1024]:
                c={x['method']:x['mse'] for x in cells if x['family']==family and x['delta']==delta and x['n']==n}
                lines.append(f"| {family} | {delta} | {n} | {c['neural_adaptive']:.8f} | {c['neural_rich']:.8f} | {c['neural_fixed']:.8f} | 3 |")
    lines+=['', '## Actual diagnostics','', '| Diagnostic | Value |','|---|---:|']
    for k in ['rotation_models','selected_corrections','max_score_gram_refinement_difference','max_target_gradient_refinement_difference','max_fitted_score_center_abs','max_source_crps_refinement_difference','max_selected_solver_relative_residual','max_selected_regularized_condition','finite_paths','source_compatible_finite_paths','max_abs_finite_path_B','max_abs_nonlinear_remainder','accepted_enrichment_directions','witness_after_to_before_median','max_truth_target_refinement_difference']:
        lines.append(f'| {k} | {diagnostics[k]} |')
    lines+=['',f'Adaptive lambda counts across the 72 rotations: `{dict(lambdas)}`; selected dimensions: `{dict(dimensions)}`.', '',
            'The weighted Sobolev metric refinement difference is 1.368362873580413 in matrix spectral norm (4096 vs 2048 nodes). A small linear solver residual is not an integral-error bound, and a finite witness is not a full-class residual bound. The near-diagonal witness plot shows that enrichment often changed its own stabilized residual indicator very little. The finite-path remainder is retained rather than assumed zero.', '',
            '## Figures and raw artifacts','',
            'All figures were rendered from the saved pilot outputs and visually inspected. Each is available as PNG and standalone SVG; the underlying data and raw-file hash are saved.', '',
            f'1. [Interaction error by family, separation and sample size](../../../{run}/evaluation/figure1_interaction_error.png) — means and every replicate; [SVG](../../../{run}/evaluation/figure1_interaction_error.svg).',
            f'2. [Same-fit ablations and CPU costs](../../../{run}/evaluation/figure2_ablation_cost.png) — shared work fully charged to each comparator; [SVG](../../../{run}/evaluation/figure2_ablation_cost.svg).',
            f'3. [Residual witnesses and finite-path nonlinearity](../../../{run}/evaluation/figure3_witness_nonlinearity.png) — lower witnesses, not uncertainty upper bounds; [SVG](../../../{run}/evaluation/figure3_witness_nonlinearity.svg).', '',
            f'[All 240 raw paired rows](../../../{run}/evaluation/raw_results.csv), [machine-readable summary](../../../{run}/evaluation/summary.json), [diagnostics](../../../{run}/evaluation/diagnostics.json), [plotted numerical data](../../../{run}/evaluation/figure_data.json), [prediction hashes](../../../{run}/evaluation/prediction_hashes.json). Per-dataset source folds, model parameters, source losses, selections, directions and endpoint traces are in `{run}/predictions/dNNN/`.', '',
            '## Failures, costs and limits','',
            'No development/final fit or direct-target prediction failed. The intentional one-second watchdog test returned 124 and its spawned child was confirmed absent. The initial discovery mistake appended historical job 27 before the new prompt prohibition was read; this was disclosed, preserved and also charged to the fresh allowance. Default filesystem image viewing failed; budgeted reads of the existing images succeeded. Scalar-logging warnings remain in validation/run records.', '',
            f"Final fitter totals (24 datasets, all methods/shared work counted once) were {diagnostics['fitter_total_cpu_s']:.6f} CPU seconds and {diagnostics['fitter_total_wall_s']:.6f} wall seconds. These omit wrapper/discovery/development/reporting costs, which the full ledger includes. Summing the table's fully charged method costs would deliberately count shared fits repeatedly.", '',
            'The implementation computes the complete 36-dimensional Gram bank before enrichment, including for the shared fixed control. Thus it has no demonstrated scalable saving from avoiding that bank; the fixed method also pays for geometry it would not strictly need in a separately optimized implementation. Neural profile misspecification, source-objective quadrature, heuristic lambda selection, restricted direction coverage and incomplete endpoint optimization limit interpretation. More seeds cannot resolve the present missing theorem or prove the biological assumptions.', '',
            'See [the proof gaps](theorem.md), [source comparison](novelty.md), [static review](review.md), [resource accounting](resources.md) and [actual commands](commands.sh).']
    with (docs/'results.md').open('x') as f:f.write('\n'.join(lines)+'\n')
    print(json.dumps(diagnostics,indent=2))


if __name__=='__main__':main()
