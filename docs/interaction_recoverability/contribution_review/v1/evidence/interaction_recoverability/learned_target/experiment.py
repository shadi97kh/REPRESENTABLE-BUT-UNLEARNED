"""Learner orchestration with source-only inputs and frozen predictions."""
import hashlib
import json
import time
from pathlib import Path
import numpy as np
from .profiles import target
from .source_fit import fit
from .geometry import metric,fit_geometry_error
from .correction import score_observations,evaluation_correction
from .enrichment import select
from .nonlinear_bias import finite_path
from .baselines import direct_target_set_search

DEFAULTS={'width':4,'learning_rate':.03,'checkpoints':[100,200,300],
          'source_nodes':64,'geometry_nodes':96,'metric_nodes':4096,
          'lambdas':[.00001,.001,.1],'enrichment_rounds':4,'search_tau':.001,
          'endpoint_steps':40,'roles':3,'selection_rule':'local displacement proxy; heuristic, not certified risk',
          'prediction_clipping':False}

def serializable(value):
    if isinstance(value,np.ndarray):return value.tolist()
    if isinstance(value,np.generic):return value.item()
    if isinstance(value,dict):return {k:serializable(v) for k,v in value.items()}
    if isinstance(value,(list,tuple)):return [serializable(v) for v in value]
    return value

def write_json(path,value):
    with Path(path).open('x') as f:json.dump(serializable(value),f,indent=2,allow_nan=False)

def stratified_roles(samples,seed):
    rng=np.random.default_rng(seed);parts=[]
    for y in samples:parts.append([x.tolist() for x in np.array_split(rng.permutation(len(y)),3)])
    return parts

def load_source(path):
    with np.load(path,allow_pickle=False) as data:
        if set(data.files)!={f'y{a}' for a in range(5)}:raise ValueError('source archive must contain only five outcome arrays')
        samples=[np.array(data[f'y{a}'],dtype=float) for a in range(5)]
    if any(y.ndim!=1 or len(y)<9 or not np.isfinite(y).all() or (y<=0).any() for y in samples):
        raise ValueError('invalid source-only outcomes')
    return samples

def fit_dataset(source_path,out,delta,seed,settings,M):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    start=time.process_time();wall=time.monotonic();samples=load_source(source_path)
    parts=stratified_roles(samples,seed+917);write_json(out/'folds.json',parts)
    predictions={k:[] for k in ['neural_shared_plugin','neural_fixed','neural_adaptive','neural_random','neural_rich','basis_shared_plugin','basis_rich']}
    fold_records=[];nuisance_cost={'neural':0.,'basis':0.}
    shared_geometry_cost={'neural':0.,'basis':0.};method_extra={}
    for fold in range(3):
        F=[y[np.array(p[fold])] for y,p in zip(samples,parts)]
        C=[y[np.array(p[(fold+1)%3])] for y,p in zip(samples,parts)]
        E=[y[np.array(p[(fold+2)%3])] for y,p in zip(samples,parts)]
        counts=[len(y) for y in E];pi=np.array(counts)/sum(counts)
        record={'fold':fold,'roles':{'fit':fold,'selection':(fold+1)%3,'evaluation':(fold+2)%3},'evaluation_counts':counts,'models':{}}
        for kind in ['neural','basis']:
            model,fitlog=fit(F,C,delta,kind,seed+fold*11+(0 if kind=='neural' else 1),settings)
            nuisance_cost[kind]+=fitlog['cpu_s'];plugin=float(target(model,delta,settings['geometry_nodes']))
            predictions[kind+'_shared_plugin'].append(plugin)
            geometry_start=time.process_time()
            geo,errors=fit_geometry_error(model,delta,pi,M,settings['geometry_nodes'])
            Cscore=score_observations(model,C,delta,geo['center'])
            shared_geometry_cost[kind]+=time.process_time()-geometry_start
            modes=['fixed','adaptive','random','rich'] if kind=='neural' else ['rich']
            selections={mode:select(geo,Cscore,counts,mode,settings,seed+100+fold) for mode in modes}
            # E enters only after all source-only choices are finalized.
            evaluation_start=time.process_time()
            Escore=score_observations(model,E,delta,geo['center'])
            shared_geometry_cost[kind]+=time.process_time()-evaluation_start
            methods={}
            for mode,solution in selections.items():
                correction,variance=evaluation_correction(Escore,solution['c'])
                prediction=plugin+correction
                if not np.isfinite(prediction):raise ArithmeticError('nonfinite corrected prediction')
                predictions[kind+'_'+mode].append(prediction)
                method_extra[kind+'_'+mode]=method_extra.get(kind+'_'+mode,0.)+solution['selection_and_construction_cpu_s']
                methods[mode]={**solution,'prediction':prediction,'correction':correction,'evaluation_variance_estimate':variance,
                               'conditional_variance_note':'per-fold estimate only; rotations are dependent'}
                if mode=='adaptive':methods[mode]['finite_path']=finite_path(model,geo,solution['c'],solution['witness_direction'],delta,C)
            record['models'][kind]={'fit':fitlog,'payload':model.payload(),'plugin':plugin,'geometry_error':errors,'methods':methods}
        write_json(out/f'fold-{fold}.json',record);fold_records.append(record)
        print(json.dumps({'output':str(out),'fold_complete':fold,'elapsed_cpu_s':time.process_time()-start}),flush=True)
    final={name:float(np.mean(values)) for name,values in predictions.items()}
    pooled={}
    for kind in ['neural','basis']:
        model,log=fit(samples,None,delta,kind,seed+80+(kind=='basis'),settings)
        pooled[kind]={'payload':model.payload(),'fit':log}
        final[kind+'_pooled_plugin']=float(target(model,delta,settings['geometry_nodes']))
        if kind=='basis':
            control=direct_target_set_search(model,samples,delta,settings)
            pooled[kind]['direct_target_set_control']=control
            final['basis_direct_set_midpoint']=control['prediction']
    method_cost={kind+'_shared_plugin':nuisance_cost[kind] for kind in ['neural','basis']}
    for name,value in method_extra.items():
        kind=name.split('_')[0];method_cost[name]=nuisance_cost[kind]+shared_geometry_cost[kind]+value
    for kind in ['neural','basis']:method_cost[kind+'_pooled_plugin']=pooled[kind]['fit']['cpu_s']
    method_cost['basis_direct_set_midpoint']=pooled['basis']['fit']['cpu_s']+pooled['basis']['direct_target_set_control']['cpu_s']
    record={'method_cpu_s':method_cost,'shared_geometry_cpu_s':shared_geometry_cost,'predictions':final,'fold_predictions':predictions,'pooled':pooled,'nuisance_cpu_s':nuisance_cost,
            'total_cpu_s':time.process_time()-start,'total_wall_s':time.monotonic()-wall,
            'source_sha256':hashlib.sha256(Path(source_path).read_bytes()).hexdigest(),
            'status':'fitted predictions frozen before evaluation-only truth access',
            'selection_access':'F and C only; E only for final correction; no evaluator imports or inputs',
            'crossfit_note':'three cyclic roles; no independence assumption across rotation estimates'}
    write_json(out/'predictions.json',record)
    return record

def prepare_metric(out,settings):
    out=Path(out)
    M=metric(settings['metric_nodes']);comparison=metric(settings['metric_nodes']//2)
    np.save(out/'metric.npy',M)
    write_json(out/'metric_accuracy.json',{'refinement_difference_spectral_norm':float(np.linalg.norm(M-comparison,2)),
        'smallest_eigenvalue':float(np.linalg.eigvalsh(M).min()),
        'status':'numerical X-metric; finite-sieve normalization, not a certified infinite-line quadrature enclosure'})
    return M
