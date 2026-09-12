"""Authorized learned-target pilot stages with exclusive outputs and frozen interfaces."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import time
import numpy as np
from .experiment import DEFAULTS,write_json,fit_dataset,prepare_metric

def load_config(run):
    path=run/'frozen_config.json';config=json.loads(path.read_text())
    for name,digest in config['learner_code_sha256'].items():
        if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest:
            raise RuntimeError('Frozen learner code changed: '+name)
    return config

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    sub=parser.add_subparsers(dest='stage',required=True)
    for stage in ['validate','development','freeze','generate','evaluate','plot','verify']:
        sub.add_parser(stage).add_argument('--run-root',required=True)
    train=sub.add_parser('train');train.add_argument('--run-root',required=True)
    train.add_argument('--start',type=int,required=True);train.add_argument('--stop',type=int,required=True)
    args=parser.parse_args();run=Path(args.run_root)
    if args.stage=='validate':
        result=subprocess.run(['python','-m','pytest','-q','tests/interaction_recoverability_learned_target'],capture_output=True,text=True)
        write_json(run/'validation.json',{'argv':['python','-m','pytest','-q','tests/interaction_recoverability_learned_target'],
                   'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        print(result.stdout);raise SystemExit(result.returncode)
    if args.stage=='development':
        from .evaluate import generate_dataset,truth_target
        dev=run/'development';dev.mkdir(exist_ok=False)
        write_json(dev/'settings.json',DEFAULTS);M=prepare_metric(dev,DEFAULTS)
        generate_dataset(dev/'source.npz',dev/'private/truth.json',family='oscillatory',delta=.03,n=1024,generator_seed=9181199,identifier='development')
        pred=fit_dataset(dev/'source.npz',dev/'fitted',.03,72001,DEFAULTS,M)
        truth=json.loads((dev/'private/truth.json').read_text())
        write_json(dev/'evaluation.json',{'target':truth['target'],'squared_errors':{k:(v-truth['target'])**2 if v is not None else None for k,v in pred['predictions'].items()},
                   'use':'development execution check only; configuration choice uses runtime, not these target errors'})
        print(json.dumps({'development_complete':True,'cpu_s':pred['total_cpu_s'],'wall_s':pred['total_wall_s']}))
    elif args.stage=='freeze':
        dev=json.loads((run/'development/fitted/predictions.json').read_text())
        ledger=[json.loads(x) for x in (run/'ledger.jsonl').read_text().splitlines()]
        remaining=7200-sum(e.get('conservative_charge_s',0) for e in ledger)-600-60
        full=24*1.5*dev['total_wall_s'];minimum=12*1.5*dev['total_wall_s']
        choice='full' if full<=remaining else 'minimum' if minimum<=remaining else 'infeasible'
        if choice=='infeasible':raise RuntimeError('Even predeclared minimum exceeds measured remaining resources')
        settings=json.loads((run/'development/settings.json').read_text())
        public=[];private=[];idx=0
        for family in ['oscillatory','localized']:
            for delta in [.1,.03]:
                for n in ([256,1024] if choice=='full' else [1024]):
                    for replicate in range(3):
                        identifier=f'd{idx:03d}'
                        public.append({'id':identifier,'delta':delta,'n':n,'initialization_seed':72100+100*idx})
                        private.append({'id':identifier,'family':family,'delta':delta,'n':n,'generator_seed':9181200+idx,'replicate':replicate})
                        idx+=1
        (run/'private').mkdir(exist_ok=False)
        write_json(run/'private/evaluation_spec.json',{'datasets':private,'profile_definitions':'evaluate.truth_profile, frozen code hash',
                   'coefficient_rule':'independent Uniform(.6,.9) for all four coefficients in each dataset'})
        code={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in Path('interaction_recoverability/learned_target').glob('*.py')}
        config={'version':1,'status':'FROZEN BEFORE FINAL DATA GENERATION AND TARGET ACCESS','plan':choice,
                'settings':settings,'datasets':public,'learner_code_sha256':code,
                'evaluation_spec_sha256':hashlib.sha256((run/'private/evaluation_spec.json').read_bytes()).hexdigest(),
                'selection_from_timing_only':{'development_wall_s':dev['total_wall_s'],'full_projection_s':full,'minimum_projection_s':minimum,'available_compute_s':remaining},
                'primary':'ordinary squared interaction error, untransformed, no clipping',
                'screening_rule':'at least 10% aggregate MSE reduction over strongest complete non-oracle control; descriptive only',
                'roles':'three stratified cyclic F,C,E rotations; source checkpoints and lambda selected without own E data',
                'methods':['neural_shared_plugin','neural_fixed','neural_adaptive','neural_random','neural_rich','basis_shared_plugin','basis_rich','neural_pooled_plugin','basis_pooled_plugin','basis_direct_set_midpoint'],
                'stops':['resource reserve','nonfinite or failed data fit retained without target-driven restart','no fits after final evaluation'],
                'truth_access':'separate evaluation module/files; no latent arrays or generator seeds passed to learner'}
        write_json(run/'frozen_config.json',config)
        destination=Path('configs/interaction_recoverability_learned_target_v1.json')
        with destination.open('x') as f:f.write((run/'frozen_config.json').read_text())
        print(json.dumps({'plan':choice,'datasets':len(public),'sha256':hashlib.sha256((run/'frozen_config.json').read_bytes()).hexdigest(),'projection_s':full if choice=='full' else minimum}))
    elif args.stage=='generate':
        config=load_config(run)
        from .evaluate import generate_final
        generate_final(run,config);print('Generated source-only files and separate evaluator truth.')
    elif args.stage=='train':
        config=load_config(run)
        if (run/'evaluation').exists():raise RuntimeError('Final targets already opened; further fits forbidden')
        if not 0<=args.start<args.stop<=len(config['datasets']):raise ValueError('invalid frozen dataset slice')
        M=np.load(run/'development/metric.npy',allow_pickle=False)
        for spec in config['datasets'][args.start:args.stop]:
            pred=fit_dataset(run/'source'/f"{spec['id']}.npz",run/'predictions'/spec['id'],spec['delta'],spec['initialization_seed'],config['settings'],M)
            print(json.dumps({'dataset_complete':spec['id'],'cpu_s':pred['total_cpu_s'],'wall_s':pred['total_wall_s']}),flush=True)
    elif args.stage=='evaluate':
        config=load_config(run)
        from .evaluate import score_final
        result=score_final(run,config);print(json.dumps({k:v for k,v in result.items() if k!='per_cell'},indent=2))
    elif args.stage=='plot':
        config=load_config(run)
        if (run/'evaluation/figure_data.json').exists():raise FileExistsError('figure outputs already exist')
        from .evaluate import plot_results
        plot_results(run);print('Three evidence figures exported as PNG and SVG.')
    elif args.stage=='verify':
        config=load_config(run);manifest=json.loads((run/'preservation.json').read_text())
        changed=[name for name,digest in manifest.items() if not Path(name).exists() or hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest]
        result={'preservation_count':len(manifest),'changed_historical_files':changed,'frozen_code_verified':True,
                'tracked_diff_empty':subprocess.run(['git','diff','--quiet']).returncode==0,
                'prediction_count':len(list((run/'predictions').glob('*/predictions.json'))),
                'figure_count':len(list((run/'evaluation').glob('figure*.svg')))}
        write_json(run/'verification.json',result);print(json.dumps(result));assert not changed

if __name__=='__main__':main()
