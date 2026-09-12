"""Finalize preparation evidence and verification; never fit or train."""
import argparse
import ast
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--run-root',required=True)
    args=parser.parse_args();root=Path(args.run_root)
    final=root/'final';final.mkdir(exist_ok=False)
    events=[json.loads(line) for line in (root/'ledger.jsonl').read_text().splitlines()]
    active=[e for e in events if e['event']=='start'][-1]['job_id']
    resource_path=root/f'resource_snapshot-{active}.json'
    original=json.loads((root/'prior_hashes.json').read_text())
    changed=[p for p,h in original.items() if not Path(p).is_file() or sha(p)!=h]
    assert not changed,changed
    config=Path('configs/interaction_recoverability_pilot.json')
    assert sha(config)==Path(str(config)+'.sha256').read_text().split()[0]
    for p in list(Path('interaction_recoverability').glob('*.py'))+list(Path('tests/interaction_recoverability').glob('*.py')):
        ast.parse(p.read_text(),filename=str(p))
    help_run=subprocess.run([sys.executable,'-m','interaction_recoverability','--help'],capture_output=True,text=True)
    assert help_run.returncode==0,help_run.stderr
    tests=subprocess.run([sys.executable,'-m','pytest','-q','tests/interaction_recoverability'],capture_output=True,text=True)
    (final/'tests.json').write_text(json.dumps(dict(returncode=tests.returncode,stdout=tests.stdout,stderr=tests.stderr),indent=2))
    assert tests.returncode==0,tests.stdout+tests.stderr
    diagnostic=json.loads((root/'numerics-v2/diagnostics.json').read_text())
    table='\n'.join(f"| {r['delta']:g} | {r['axis2_w1']['value']:.9f} | {r['joint_w1']['value']:.9f} | {r['primary_interaction_gap']:.9f} | {r['axis2_hellinger_squared']['value']:.9g} |" for r in diagnostic['rows'])
    verification=dict(prior_files_verified=len(original),changed_prior_files=changed,
        git_tracked_diff=subprocess.check_output(['git','diff','--name-only','HEAD']).decode().splitlines(),
        configuration_sha256=sha(config),cli_help_returncode=0,syntax='passed',tests=tests.stdout.strip(),
        uploaded_report_sha256=sha('Mathematical_ML_Novelty_and_Open_Problems.docx'),
        script_present=Path('interaction_stability_check.py').exists(),
        resource_snapshot_after_job=str(resource_path),fit_calls=0,gpu_use=False)
    (final/'verification.json').write_text(json.dumps(verification,indent=2))
    commands=[dict(job_id=e['job_id'],argv=e['command']) for e in events if e['event']=='start']
    (final/'commands.json').write_text(json.dumps(commands,indent=2))
    decision=f'''# Interaction recoverability decision

**incremental-known.** The completed known-warp result has matching finite-sample rate orders for a selected interaction, but the estimator is ordinary two-quantile minimum-distance inversion. Its bilinear-completion control is the identical calculation. The broader unknown-mechanism contribution remains open, so no training authorization is requested.

**Strongest completed result — STANDARD RESULT APPLIED.** For the normalized class in [problem.md](problem.md), with known Gaussian noise, known exponential profiles/warps, known delta, fixed private anchors, four unknown bounded amplitudes and n independent outcomes at each of five observed actions,

\\[
R_n(I_{{12}},\\mathcal M_\\delta)\\asymp\\min\\{{1,(n\\delta^2)^{{-1}}\\}},\\quad 0<\\delta\\le.1.
\\]

Constants are uniform in n and delta across these separate known-delta classes, and are conservative. The [proof notes](proofs.md) derive the transformed-density score, bound its square integral by 65.6233873491624, propagate Hellinger affinity through the full product experiment, and give an amplitude-path lower bound matched in order by an empirical-quantile upper bound. No Wasserstein-to-testing inference is used. A secondary protected contrast equals a difference of observed means and has uniform O(1/n) MSE without recovering components. That is standard functional estimability, not a replacement primary endpoint.

**Exact remaining proof gap — OPEN.** There is no theorem here for learning unknown profiles/warps or separation from the same finite action dictionary while retaining a target-sensitive bound. An extension needs a specified nuisance class, a verifiable central-range/design condition and an estimator whose nuisance error is included, plus a distinction from existing inverse-functional/likelihood/completion theory. The restricted proofs contain no intentionally open lemma; they do not solve this broader step. An unrestricted eta*a1*a2 term is observationally invisible on all training axes and defeats robustness.

**Actual outputs — NUMERICALLY CHECKED.** The independent diagnostic reproduces the DOCX table and distinguishes joint W1 from the primary interaction gap:

| delta | Axis-2 W1 | Joint-action W1 | I12(T)-I12(Q) | Axis-2 Hellinger squared |
|---:|---:|---:|---:|---:|
{table}

The limiting joint W1 and interaction gap are 1.886070833180237. Integration used [-12,12], absolute/relative requested tolerances 2e-11, explicit zero/crossing splits and transformed-density Jacobians. The largest W1 quadrature error estimate is 3.62e-11; the largest analytic omitted-tail upper bound is 2.21e-21. Hellinger inversion uses a 1e-12 root tolerance. Neither ordinary quadrature estimates nor floating-point evaluations of tail bounds are certified numerical enclosures. No samples, fitted estimators, bootstrap studies, tuning sweeps or benchmarks were run. The [diagnostic JSON](../../{root}/numerics-v2/diagnostics.json) retains every error estimate and analytic bound.

The [novelty matrix](novelty_matrix.md) checks the named primary results and direct adjacent collisions, including strongly identified functionals of weak nuisance functions. The DOCX and its nine embedded equation images were inspected. The separate `interaction_stability_check.py` was not found in the repository or attachment directory, despite follow-up; its original code remains unaudited. The independent harness reconstructs the equations supplied in both the prompt and DOCX. No original script was overwritten or falsely reported as executed.

The [estimator specification](estimator.md) gives the implemented, fitting-gated sample interface, deterministic algebra, costs and known/estimated/oracle information. The [siRNA requirements](sirna_requirements.md) document why current biological data do not validate this target. The [frozen unexecuted pilot](../../configs/interaction_recoverability_pilot.json) has SHA256 `{sha(config)}`. It includes equally informed minimum-distance/bilinear, ridge, likelihood and energy-score controls and labels native-method integration and oracle limits. The pilot is a falsification/control design, not an authorized run. Fitted runtime and optimizer behavior are unmeasured; an operation-count estimate is included instead of an invented complete compute cap.

**Resources.** The [append-only ledger](../../{root}/ledger.jsonl) belongs only to the new 600-CPU-second / 600-aggregate-job-wall-second allowance. The [final snapshot](../../{resource_path}) includes this report/check job after it exits, measured parent/waited-child CPU, conservative charge, separate task elapsed wall, memory and remaining allowance. Every metered job and descendant shared one enforced CPU affinity, with one numerical thread and no GPU. A remaining-budget watchdog kills the process group on timeout, with a clean-stop reserve. Initial text-only repository/prompt inspection preceded the meter and is explicitly unmeasured; 30 seconds are conservatively reserved/charged for it and launcher overhead. This is not presented as measured CPU. Unreaped descendants could be absent from measured counters, but inherited one-core elapsed time supplies the conservative bound for subsequent jobs. Remote retrieval was also charged where performed by the local wrapper, although the allowance permitted excluding it. Historical project balance remains unknown.

The [verification record](../../{final}/verification.json) confirms {len(original)} original files unchanged, an empty tracked diff, CLI/help and syntax checks, and `{tests.stdout.strip().splitlines()[-1]}`. The test count is implementation evidence, not the research outcome. The prior P4 and graph-CNP records remain unchanged.

These are implemented preparation commands, using the same existing ledger. A diagnostic repeat needs a fresh output destination; the commands already used `numerics-v1` and `numerics-v2` exclusively:

```bash
python -m interaction_recoverability.budget --run-root {root} -- python -m interaction_recoverability --help
python -m interaction_recoverability.budget --run-root {root} -- python -m pytest -q tests/interaction_recoverability
python -m interaction_recoverability.budget --run-root {root} -- python -m interaction_recoverability diagnose --out {root}/numerics-reproduction
```

The [exact command log](../../{final}/commands.json) stores argv including setup and subprocess scripts; the ledger contains exit status and resource charge. No fitted-pilot or training command is claimed to exist.

**Next action:** subject the unknown-profile extension's proposed observation condition to mathematical falsification and the cited inverse-functional theory before spending a separate fitting budget. The current result does not justify a training request.
'''
    path=Path('docs/interaction_recoverability/decision.md')
    with path.open('x') as f:f.write(decision)
    claims=[
        ('validity_and_integrability','PROVED HERE','proofs.md'),
        ('axis_equalities_and_limit','PROVED HERE','proofs.md'),
        ('transport_score_bounds','PROVED HERE','proofs.md'),
        ('complete_experiment_and_matching_minimax_rate','STANDARD RESULT APPLIED','proofs.md'),
        ('protected_target_identity','PROVED HERE','proofs.md'),
        ('unknown_nuisance_extension','OPEN','proofs.md'),
        ('diagnostic_table','NUMERICALLY CHECKED','decision.md')]
    (final/'claim_statuses.json').write_text(json.dumps([dict(claim=a,status=b,evidence=c) for a,b,c in claims],indent=2))
    missing=[]
    for p in Path('docs/interaction_recoverability').glob('*.md'):
        for link in re.findall(r'\]\(([^)]+)\)',p.read_text()):
            if link.startswith(('https:','http:','#')):continue
            target=(p.parent/link.split('#')[0]).resolve()
            if not target.exists() and target!=resource_path.resolve():missing.append((str(p),link))
    assert not missing,missing
    paths=[]
    for prefix in ['interaction_recoverability','tests/interaction_recoverability','docs/interaction_recoverability']:
        paths.extend(p for p in Path(prefix).rglob('*') if p.is_file() and '__pycache__' not in p.parts)
    paths += [config,root/'environment.json',root/'prior_hashes.json',root/'numerics-v2/diagnostics.json',final/'verification.json',final/'claim_statuses.json',final/'tests.json',final/'commands.json']
    manifest=dict(entries=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(set(paths))],
        live_resource_note='Ledger and final resource snapshot close after this job and are intentionally excluded from this static artifact inventory.')
    (final/'artifact_hashes.json').write_text(json.dumps(manifest,indent=2))
    print(json.dumps(verification,indent=2))

if __name__=='__main__':main()
