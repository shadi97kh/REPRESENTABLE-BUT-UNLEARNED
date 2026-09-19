"""Finalize existing-artifact provenance and command journal; no numerical experiments."""
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
from pathlib import Path
import shlex


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--run-root',required=True);p.add_argument('phase',choices=['report'])
    args=p.parse_args();run=Path(args.run_root);docs=Path(__file__).resolve().parent
    if (run/'final_checks.json').exists():raise FileExistsError('Finalization already completed')
    config=json.loads((run/'frozen_config.json').read_text())
    assert all(sha(path)==digest for path,digest in config['learner_code_sha256'].items())
    assert sha(run/'private/evaluation_spec.json')==config['evaluation_spec_sha256']
    saved=json.loads((run/'evaluation/prediction_hashes.json').read_text())
    assert len(saved)==24
    for name,digest in saved.items():
        prediction=run/'predictions'/name/'predictions.json';assert sha(prediction)==digest
        obj=json.loads(prediction.read_text());assert sha(run/'source'/f'{name}.npz')==obj['source_sha256']
    original=json.loads((run/'preservation.json').read_text())
    changed=[name for name,digest in original.items() if not Path(name).exists() or sha(name)!=digest]
    assert not changed
    events=[json.loads(line) for line in (run/'ledger.jsonl').read_text().splitlines()]
    starts=[e for e in events if e['event']=='start'];ends={e['job_id']:e for e in events if e['event']=='end'}
    current=[s for s in starts if s['job_id'] not in ends];assert len(current)==1
    records=[];stage_lines=['#!/usr/bin/env bash','set -eu','',
        '# EXECUTED command journal, through completed job '+str(max(ends)),
        '# The existing outputs are exclusive and intentionally refuse replay.',
        '# This file contains the tested operational forms, not authorization to rerun fits.',
        '# All discovery/read/image/report argv and statuses: command_records.json; ledger.jsonl is authoritative.',
        '# Wrapper prefixes below are reconstructed from recorded phase/timeout; child argv are exact.',
        '', '# Initialization ran successfully before job 1:',
        shlex.join(['python','-m','interaction_recoverability.learned_target.budget','--run-root',str(run),'--init']), '']
    for s in starts:
        if s['job_id'] not in ends:continue
        e=ends[s['job_id']];wrapper=['python','-m','interaction_recoverability.learned_target.budget','--run-root',str(run),'--phase',s['phase'],'--job-seconds',str(s['timeout_s']),'--']
        rc=124 if e['timed_out'] else e['returncode']
        rec={'job_id':s['job_id'],'exact_child_argv':s['command'],'reconstructed_wrapper_argv':wrapper+s['command'],
             'wrapper_exit_code':rc,'child_returncode':e['returncode'],'timed_out':e['timed_out'],
             'start_utc':s['utc'],'end_utc':e['utc'],'resource_snapshot':str(run/f"resource_snapshot-{s['job_id']}.json")}
        records.append(rec)
        operational=('interaction_recoverability.learned_target' in s['command'] or 'pytest' in s['command'] or e['timed_out'] or any('compile_report.py' in x for x in s['command']))
        if operational:
            stage_lines.extend([f"# Job {s['job_id']}: EXECUTED; exit {rc}; snapshot resource_snapshot-{s['job_id']}.json",shlex.join(wrapper+s['command']),''])
    stage_lines+=['# Report finalization itself runs after this journal cutoff; its exact command and',
                  '# eventual exit status are in ledger.jsonl and the final resource snapshot.',
                  '# No future optional fits are listed or authorized here.']
    with (docs/'commands.sh').open('x') as f:f.write('\n'.join(stage_lines)+'\n')
    with (docs/'command_records.json').open('x') as f:json.dump(records,f,indent=2)
    # A debit, never a reset/refund: cover text/file/image tool overhead outside shell timing.
    extra={'utc':datetime.now(timezone.utc).isoformat(),'event':'overhead','conservative_charge_s':120.,
           'measured_cpu_s':None,'note':'Additional conservative debit for unmetered implementation/document-authoring and image-tool/file I/O overhead, including final resource-note creation. Distinct from measured CPU/job wall; never increases allowance.'}
    with (run/'ledger.jsonl').open('a') as f:f.write(json.dumps(extra)+'\n');f.flush();os.fsync(f.fileno())
    results=run/'evaluation';checks={'allowance_start_utc':events[0]['utc'],'status':'all existing-artifact checks passed','protected_files':len(original),'changed_protected_files':changed,
        'frozen_code_unchanged':True,'private_spec_commitment_matches':True,'prediction_and_source_hashes_match':True,
        'datasets':len(saved),'png_figures':len(list(results.glob('figure*.png'))),'svg_figures':len(list(results.glob('figure*.svg'))),
        'watchdog':json.loads((run/'watchdog_check.json').read_text()),'historical_ledger_exception':'Discovery job27, disclosed and preserved; no subsequent historical edits',
        'last_running_job':current[0]['job_id'],'final_snapshot':str(run/f"resource_snapshot-{current[0]['job_id']}.json"),
        'command_journal_cutoff_job':max(ends),'extra_unmetered_overhead_debit_s':120.,
        'static_review':'docs/interaction_recoverability/learned_target_v1/review.md; executing-agent self-review, no independent-agent claim'}
    with (run/'final_checks.json').open('x') as f:json.dump(checks,f,indent=2)
    # Manifest intentionally excludes live ledger/snapshots and the final resource note,
    # which is written from the wrapper's completed snapshot after this process returns.
    files=list(docs.glob('*'))+list(results.glob('*'))+[run/'frozen_config.json',run/'final_checks.json']
    manifest={str(x):sha(x) for x in files if x.is_file()}
    with (run/'artifact_hashes.json').open('x') as f:json.dump(manifest,f,indent=2)
    print(json.dumps(checks,indent=2))


if __name__=='__main__':main()
