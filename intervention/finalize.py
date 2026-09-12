"""Package the completed preparation; never starts training or network retrieval."""
from pathlib import Path
import importlib.metadata
import json
import platform
import statistics
import subprocess
from .accounting import sha, write_json

ROOT = Path(__file__).resolve().parents[1]


def main():
    regions = []
    for folder in ('data/intervention_sources', 'runs/intervention'):
        for p in sorted((ROOT/folder).rglob('*.json')):
            d = json.loads(p.read_text())
            if isinstance(d, dict) and isinstance(d.get('resources'), dict) and 'cpu_core_hours' in d['resources']:
                regions.append(dict(path=str(p.relative_to(ROOT)), **d['resources']))
    profile = json.loads((ROOT/'runs/intervention/preparation-v2/profile.json').read_text())
    online = json.loads((ROOT/'runs/intervention/online_cost_diagnostic-v1.json').read_text())
    medians = {
        name: statistics.median(m[name]['process_cpu_s'] + m[name]['waited_child_cpu_s'] for m in online['measurements'])
        for name in ('direct_uncached_eight_states','direct_warm_cache','untrained_operator_inference')
    }
    summary = dict(
        date='2026-09-12', scope='New intervention preparation only; measured regions, not a complete session ledger',
        measured_regions=regions,
        measured_cpu_core_hours=sum(r['cpu_core_hours'] for r in regions),
        measured_process_cpu_s=sum(r['process_cpu_s'] for r in regions),
        measured_waited_child_cpu_s=sum(r['waited_child_cpu_s'] for r in regions),
        sum_instrumented_region_wall_s=sum(r['wall_s'] for r in regions),
        wall_sum_limit='This sum is not end-to-end elapsed time; orchestration gaps and uninstrumented requests/editing are excluded.',
        maximum_reported_parent_peak_rss_kib=max(r['peak_rss_kib'] for r in regions),
        memory_limit='Parent process high-water marks, not a simultaneous process-tree maximum; child RSS not included.',
        gpu_hours=0, optimizer_steps=0, predictor_training_steps=0,
        explainer_training_steps=0, held_out_experiments=0,
        historical_remaining_allowance=None,
        prospective_profile=profile['estimates'],
        online_cpu_medians_s=medians,
        observed_online_break_even=None,
        break_even_limit='Fixed analytic correctness fixture and untrained forward pass only; proposed online CPU exceeds exact evaluation even before transcript queries/training.',
        unmeasured=['tool/editor and other shell CPU outside recorded regions',
                    'unwrapped source retrieval, parsing and orchestration gaps',
                    'source coefficient materialization separately from fixture graph/query timing',
                    'Adam steps and optimizer memory; actual batching and auxiliary-loss costs',
                    'main teacher cache serialization/I/O; baseline selection and tree/SPEX fitting',
                    'published delivery inference and all model-specific descriptors/3D construction',
                    'predictor fitting, training convergence and any GPU variant'],
        versions={name:importlib.metadata.version(name) for name in ('numpy','torch','scipy','pandas','rdkit','scikit-learn','pytest')},
        python=platform.python_version(), platform=platform.platform(),
    )
    write_json(ROOT/'runs/intervention/resource_summary-v1.json', summary)

    text = f'''# Measured resources and prospective accounting

No predictor or explainer training, model zoo, policy optimization or main teacher-cache generation ran. Correctness fixtures, metadata audits and a few forward/backward checks ran on CPU with no optimizer updates. GPU use is zero. The historical balance under the old cap remains unknown and has not been reset.

The [resource summary](../../runs/intervention/resource_summary-v1.json) records {len(regions)} distinct instrumented regions, including earlier audit/profile versions and the failed intermediate test invocation. Their measured CPU total is **{summary['measured_cpu_core_hours']:.6f} core-hours**, computed from **{summary['measured_process_cpu_s']:.3f} parent process CPU seconds plus {summary['measured_waited_child_cpu_s']:.3f} waited-child CPU seconds**. Nested component profiles are excluded from this sum to avoid double counting. This is an incomplete session ledger: editor/tool work, unwrapped commands, some retrieval and orchestration are unmeasured. Their unknown cost is not zero, and this subtotal cannot establish a remaining historical allowance.

The sum of instrumented region durations is {summary['sum_instrumented_region_wall_s']:.2f} seconds; it is not end-to-end elapsed time. Maximum reported parent RSS is {summary['maximum_reported_parent_peak_rss_kib']/1024:.1f} MiB; this is not simultaneous process-tree memory. The final operator preparation reached {profile['resources']['peak_rss_kib']/1024:.1f} MiB. CPU numerical checks fixed OMP/BLAS/PyTorch threads to one and used no workers; retrieval thread settings are retained individually and are not all one. Process CPU includes all parent threads. Waited children, including PDF extraction and tests, are counted through `RUSAGE_CHILDREN`.

The final [profile](../../runs/intervention/preparation-v2/profile.json) measured 80 exact analytic fixture queries with graph construction, five untrained operator forward/backward episodes with unchanged weights, five small molecules and five actual AGILE molecule reconstructions. The operator has 17,409 parameters. First-call startup is separated from the four warm repeats, then charged once per future run. The earlier profile is preserved, but v2 is the estimate used here.

| Prospective component in frozen CPU screen | Measured basis | Estimate or missing measurement |
|---|---|---|
| Source-predictor training | Analytic coefficient functions | Zero optimizer steps by construction; coefficient creation is not separately timed |
| Complete analytic teacher tables | 21,504 queries: 16,384 train, 1,024 development, 4,096 test | {profile['estimates']['query_cache_cpu_hours_linear']:.6f} CPU core-hours by linear scaling from the 80-query fixture, including tiny graph construction; cache I/O unmeasured |
| Three loss arms × three initializations | 13,500 future optimizer steps × eight episodes = 108,000 forward/backward episodes | Warm median projection {profile['estimates']['core_hours_from_measured_median']:.3f} CPU core-hours; independently measured wall projection {profile['estimates']['wall_hours_from_measured_median']:.3f} hours |
| Planning range for those forward/backward episodes | Warm minimum to 3× warm maximum | {profile['estimates']['planning_core_hour_range'][0]:.3f}–{profile['estimates']['planning_core_hour_range'][1]:.3f} core-hours; heuristic range, not a bound or total budget |
| Cold operator startup | First observed CPU call | {profile['estimates']['cold_start_cpu_s']:.6f} CPU seconds per run, additional to the warm projection |
| Five actual lipid feature rebuilds | Eight atom features, two stereo channels, four descriptors | {profile['profiles']['five_actual_lipid_feature_rebuilds']['process_cpu_s']:.6f} CPU seconds total; not Mordred, fingerprints or TransMA 3D timing |
| Adam, actual eight-episode execution, auxiliary losses | No optimizer run | Update runtime/state memory unmeasured; heuristic multiplier does not replace measurement |
| Baseline fitting and development selection | Interfaces only, plus tiny correctness solves | Full tree/GP/LASSO/SPEX tuning and integration runtime unmeasured |
| Delivery prediction and source-model fitting | Two checkpoints inspected for tensor metadata only | Full construction, prediction, fitted membership and training costs unresolved |
| GPU | CPU-only configuration | Zero in this proposed screen; a GPU variant needs a fresh profile and prospective allowance |

The future schedule is only a concrete estimate: it was not executed and is not a requested or approved compute budget. Its optimizer, integration and data limitations prevent presenting the partial projection as a complete cap. No training command is supplied. The scientific decision is incremental known method, so preparation does not recommend spending compute simply to search for a favorable result.

The [online diagnostic](../../runs/intervention/online_cost_diagnostic-v1.json) uses one fixed analytic fixture and five timing repetitions. Median CPU times are **{1000*medians['direct_uncached_eight_states']:.3f} ms** for all eight uncached exact states, **{1000*medians['direct_warm_cache']:.3f} ms** for their warm oracle-cache lookup path and **{1000*medians['untrained_operator_inference']:.3f} ms** for an untrained operator forward pass with a transcript already available. Shared graph construction is excluded from all three. This is a timing diagnostic, not a held-out accuracy result or trained amortization comparison. Direct response-table array lookup could be cheaper still.

For `extra_upfront / (baseline_online - proposed_online)`, the measured denominator is nonpositive here, so there is no finite efficiency break-even in this observed regime. Adding transcript acquisition and training cannot repair that denominator. Stop this tiny-task efficiency claim. The measurement does not determine costs for an expensive published predictor, the full 1,200-lipid library, or a larger intervention space; those remain unprofiled.

The accounting plan for any separately authorized execution is to record stage-start/end monotonic time and parent/child CPU, set and log numerical threads and workers, wait for all child processes before taking final counters, and monitor a process tree when workers persist. Record peak process-tree RSS separately from parent high-water marks. Charge every teacher query, cache build/read/write, complete feature rebuild, source-model fit, training step, baseline fit and inference. GPU use, if introduced, requires device identity, number of devices and explicitly charged occupied device time. Enforce a newly specified prospective cap with stage checks and an external watchdog; none is inferred from the old ledger.

For the primary curve, add method-specific setup/source generation/teacher/training/preprocessing cost amortized over 1,000 deployments to actual output-completion times. Each method pays for its own required work; report total project spend as well as method allocation, and disclose any shared cached work instead of making it free. Cold exact caches pay to populate; warm reuse is a separate sensitivity. Report CPU core-hours and wall-time curves separately because parallel wall savings need not save CPU. Retain zero predictions before the first output and the latest completed output at timeouts. The frozen endpoint and stop rules are in [configuration v2](../../configs/intervention_screen_v2.json).
'''
    with (ROOT/'docs/intervention/resources.md').open('x') as f: f.write(text)

    checks=[]
    for folder in ('data/intervention_sources','runs/intervention','docs/intervention'):
        for sidecar in sorted((ROOT/folder).rglob('*.sha256')):
            p=Path(str(sidecar)[:-7])
            expected=sidecar.read_text().split()[0]
            checks.append(dict(path=str(p.relative_to(ROOT)),matches=sha(p)==expected))
    for p in sorted((ROOT/'configs').glob('intervention_screen_v*.json')):
        checks.append(dict(path=str(p.relative_to(ROOT)),matches=sha(p)==Path(str(p)+'.sha256').read_text().split()[0]))
    legacy=[]
    cycle=json.loads((ROOT/'runs/revision/p4_gpu-20260911-171200-784773-cycle.json').read_text())
    for item in cycle['artifacts'].values():
        legacy.append(dict(path=item['path'],matches=sha(ROOT/item['path'])==item['sha256']))
    for name in ('docs/p4_actual_model_executed-20260911-171200-784773-final.md',
                 'docs/p4_actual_model_executed-20260911-171200-784773.md'):
        legacy.append(dict(path=name,matches=sha(ROOT/name)==(ROOT/(name+'.sha256')).read_text().split()[0]))
    tracked=subprocess.run(['git','diff','--name-only','HEAD'],cwd=ROOT,capture_output=True,text=True,check=True).stdout.splitlines()
    check=dict(new_sidecars=checks,legacy_p4_artifacts=legacy,
               all_new_sidecars_match=all(c['matches'] for c in checks),
               all_legacy_p4_hashes_match=all(c['matches'] for c in legacy),
               tracked_diff_paths=tracked,
               evidence_file_sha256=sha(ROOT/'GNN_RNA_Data_Research_Evidence.json'),
               tests='runs/intervention/tests-v4.json',
               final_configuration='configs/intervention_screen_v2.json',
               final_configuration_sha256=sha(ROOT/'configs/intervention_screen_v2.json'))
    assert check['all_new_sidecars_match'] and check['all_legacy_p4_hashes_match']
    write_json(ROOT/'runs/intervention/integrity-v1.json',check)
    paths=[]
    for folder in ('intervention','tests/intervention','docs/intervention','runs/intervention','data/intervention_sources'):
        paths.extend(p for p in (ROOT/folder).rglob('*') if p.is_file() and '__pycache__' not in p.parts and not p.name.endswith('.sha256'))
    paths.extend((ROOT/'configs').glob('intervention_screen_v*.json'))
    write_json(ROOT/'runs/intervention/bundle-v1.json',dict(
        status='PREPARATION_COMPLETE; incremental_known_method',
        training_steps=0,gpu_hours=0,held_out_evidence=False,
        final_configuration=check['final_configuration'],
        entries=[dict(path=str(p.relative_to(ROOT)),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(set(paths))],
        inventory_excludes='its own file and sidecar; Python bytecode',
    ))
    print(json.dumps(dict(measured_cpu_core_hours=summary['measured_cpu_core_hours'],
                         regions=len(regions),new_sidecars=len(checks),legacy_p4=len(legacy),
                         tracked_diff_paths=tracked,configuration_sha256=check['final_configuration_sha256']),indent=2))


if __name__=='__main__': main()
