"""Evaluation-only truth generation, post-freeze scoring, and evidence figures.

Learner modules never import this module or receive these simulator objects/seeds.
"""
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import numpy as np
from scipy.special import roots_legendre

def dump(path,value):
    with Path(path).open('x') as f:json.dump(value,f,indent=2,allow_nan=False)

def truth_profile(x,family,j):
    if family=='oscillatory':
        amplitude,frequency=[(.045,1.2),(-.035,1.7)][j]
        mod=np.sin(frequency*x)
    else:
        amplitude,width,center=[(.045,.8,.5),(-.04,1.,-.75)][j]
        mod=np.exp(-(width*(x-center))**2)-np.exp(-(width*center)**2)
    return np.exp(x)*(1+amplitude*mod)

def truth_response(z,a,delta,coef,family):
    hz=z+delta*z*np.sqrt(1+z*z)
    return truth_profile(z+coef[0,0]*a[0]+coef[1,0]*a[1]+a[2],family,0)+truth_profile(hz+coef[0,1]*a[0]+coef[1,1]*a[1]+a[3],family,1)

def truth_target(coef,delta,family,nodes=256):
    x,w=roots_legendre(nodes);z=12*x;w=12*w*np.exp(-z*z/2)/math.sqrt(2*math.pi)
    actions=[([1,1,0,0],1),([1,0,0,0],-1),([0,1,0,0],-1),([0,0,0,0],1)]
    return float(sum(sign*w@truth_response(z,a,delta,coef,family) for a,sign in actions))

def generate_dataset(source_path,truth_path,*,family,delta,n,generator_seed,identifier):
    source_path=Path(source_path);truth_path=Path(truth_path)
    if source_path.exists() or truth_path.exists():raise FileExistsError('dataset outputs already exist')
    rng=np.random.default_rng(generator_seed);coef=rng.uniform(.6,.9,size=(2,2))
    actions=np.vstack([np.zeros(4),np.eye(4)])
    outcomes={f'y{i}':truth_response(rng.normal(size=n),a,delta,coef,family) for i,a in enumerate(actions)}
    source_path.parent.mkdir(parents=True,exist_ok=True);truth_path.parent.mkdir(parents=True,exist_ok=True)
    os.chmod(truth_path.parent,0o700)
    np.savez(source_path,**outcomes)
    I=truth_target(coef,delta,family);I2=truth_target(coef,delta,family,512)
    dump(truth_path,{'id':identifier,'family':family,'delta':delta,'n':n,'generator_seed':generator_seed,
                     'coefficients':coef.tolist(),'target':I2,'target_quadrature_difference':abs(I-I2),
                     'latent_observations_stored':False,'source_sha256':hashlib.sha256(source_path.read_bytes()).hexdigest()})
    os.chmod(truth_path,0o600)

def generate_final(run,config):
    run=Path(run);private=json.loads((run/'private/evaluation_spec.json').read_text())
    for spec in private['datasets']:
        generate_dataset(run/'source'/f"{spec['id']}.npz",run/'private'/f"{spec['id']}.json",family=spec['family'],delta=spec['delta'],n=spec['n'],generator_seed=spec['generator_seed'],identifier=spec['id'])
    dump(run/'generation_complete.json',{'datasets':len(private['datasets']),'config_sha256':hashlib.sha256((run/'frozen_config.json').read_bytes()).hexdigest()})

def score_final(run,config):
    run=Path(run);out=run/'evaluation';out.mkdir(exist_ok=False)
    rows=[];prediction_hashes={};failed=[]
    for spec in config['datasets']:
        identifier=spec['id'];p=run/'predictions'/identifier/'predictions.json'
        if not p.exists():
            failed.append(identifier);continue
        prediction_hashes[identifier]=hashlib.sha256(p.read_bytes()).hexdigest()
    # Freeze all prediction hashes before reading any final target.
    dump(out/'prediction_hashes.json',prediction_hashes)
    for spec in config['datasets']:
        identifier=spec['id']
        if identifier in failed:continue
        pred=json.loads((run/'predictions'/identifier/'predictions.json').read_text())
        truth=json.loads((run/'private'/f'{identifier}.json').read_text())
        if truth['source_sha256']!=pred['source_sha256']:raise RuntimeError('source provenance mismatch')
        for method,value in pred['predictions'].items():
            rows.append({'id':identifier,'family':truth['family'],'delta':truth['delta'],'n':truth['n'],
                         'method':method,'target':truth['target'],'prediction':value,
                         'squared_error':(value-truth['target'])**2 if value is not None else None,
                         'cpu_s_dataset_shared':pred['total_cpu_s'],'method_cpu_s_fully_charged':pred.get('method_cpu_s',{}).get(method),
                         'status':'ok' if value is not None else 'failed_direct_set_search'})
    with (out/'raw_results.csv').open('x',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    methods=sorted({r['method'] for r in rows});aggregates={}
    for method in methods:
        selected=[r for r in rows if r['method']==method and r['squared_error'] is not None]
        aggregates[method]={'mse':float(np.mean([r['squared_error'] for r in selected])) if selected else None,
                            'successful_datasets':len(selected),'failed_datasets':len(config['datasets'])-len(selected)}
    complete_controls=[m for m in methods if m!='neural_adaptive' and aggregates[m]['failed_datasets']==0]
    strongest=min(complete_controls,key=lambda m:aggregates[m]['mse']) if complete_controls else None
    advantage=1-aggregates['neural_adaptive']['mse']/aggregates[strongest]['mse'] if strongest else None
    by_cell=[]
    for family in ['oscillatory','localized']:
        for delta in [.1,.03]:
            for n in sorted({r['n'] for r in rows}):
                cell=[r for r in rows if r['family']==family and r['delta']==delta and r['n']==n]
                for method in methods:
                    vals=[r['squared_error'] for r in cell if r['method']==method and r['squared_error'] is not None]
                    if vals:by_cell.append({'family':family,'delta':delta,'n':n,'method':method,'mse':float(np.mean(vals)),'replicates':len(vals)})
    decision='empirical candidate; theory unresolved' if advantage is not None and advantage>=.1 else 'incremental-known'
    summary={'decision':decision,'aggregate':aggregates,'strongest_complete_control':strongest,
             'relative_mse_reduction_vs_strongest':advantage,'screening_threshold':.1,'per_cell':by_cell,
             'failed_datasets':failed,'datasets':len(config['datasets']),
             'limitations':['pilot evidence only, three replicates per cell','no original-class confidence certificate',
                            'no proof of sharp attainable learned-nuisance rate','finite-sieve enrichment is standard residual-driven inverse regularization']}
    dump(out/'summary.json',summary)
    return summary

def plot_results(run):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    run=Path(run);out=run/'evaluation';rows=list(csv.DictReader((out/'raw_results.csv').open()))
    rows=[{**r,'delta':float(r['delta']),'n':int(r['n']),'squared_error':float(r['squared_error']) if r['squared_error'] else None} for r in rows]
    methods=sorted({r['method'] for r in rows});colors={m:plt.get_cmap('tab10')(i) for i,m in enumerate(methods)}
    fig,axes=plt.subplots(2,2,figsize=(12,8),sharey=True)
    for ax,(family,delta) in zip(axes.flat,[(f,d) for f in ['oscillatory','localized'] for d in [.1,.03]]):
        sizes=sorted({r['n'] for r in rows});xs=np.arange(len(sizes))
        for mi,method in enumerate(methods):
            means=[]
            for x,n in zip(xs,sizes):
                vals=[r['squared_error'] for r in rows if r['family']==family and r['delta']==delta and r['n']==n and r['method']==method and r['squared_error'] is not None]
                means.append(np.mean(vals) if vals else np.nan)
                ax.scatter(x+(mi-(len(methods)-1)/2)*.035+np.linspace(-.006,.006,len(vals)),vals,color=colors[method],s=12,alpha=.45)
            ax.plot(xs+(mi-(len(methods)-1)/2)*.035,means,'o-',color=colors[method],label=method,ms=4,lw=1)
        ax.set_title(f'{family}; delta={delta}');ax.set_xticks(xs,sizes);ax.set_xlabel('Outcomes per source stratum');ax.set_yscale('log');ax.grid(alpha=.2)
    axes[0,0].set_ylabel('Squared interaction error (log scale)');axes[1,0].set_ylabel('Squared interaction error (log scale)')
    handles,labels=axes[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='lower center',ncol=3,fontsize=8)
    fig.suptitle('Held-out interaction estimation: means and every source-data replicate');fig.tight_layout(rect=[0,.14,1,.96])
    for ext in ['png','svg']:fig.savefig(out/f'figure1_interaction_error.{ext}',dpi=160)
    plt.close(fig)
    modes=['neural_shared_plugin','neural_fixed','neural_adaptive','neural_random','neural_rich']
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    for i,method in enumerate(modes):
        vals=[r['squared_error'] for r in rows if r['method']==method and r['squared_error'] is not None]
        axes[0].bar(i,np.mean(vals),color=colors[method],alpha=.6)
        axes[0].scatter(i+np.linspace(-.18,.18,len(vals)),vals,color=colors[method],s=15)
        costs=[float(r['method_cpu_s_fully_charged']) for r in rows if r['method']==method and r['method_cpu_s_fully_charged']]
        if costs:axes[1].bar(i,np.mean(costs),color=colors[method],alpha=.7)
    for ax in axes:ax.set_xticks(range(len(modes)),[m.replace('neural_','').replace('_','\n') for m in modes]);ax.grid(axis='y',alpha=.2)
    axes[0].set_yscale('log');axes[0].set_ylabel('Squared interaction error');axes[1].set_ylabel('CPU seconds, shared work charged fully per method')
    fig.suptitle('Mechanism ablation with the same nuisance fits');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'figure2_ablation_cost.{ext}',dpi=160)
    plt.close(fig)
    pairs=[];paths=[]
    for p in sorted((run/'predictions').glob('*/fold-*.json')):
        record=json.loads(p.read_text());adaptive=record['models']['neural']['methods']['adaptive']
        for tr in adaptive['trajectory']:
            if tr['accepted']:pairs.append([tr['round'],tr['witness'],tr['same_direction_residual_after']])
        for row in adaptive['finite_path']['rows']:paths.append([row['linear_prediction'],row['B'],row['nonlinear_remainder']])
    fig,axes=plt.subplots(1,2,figsize=(11,4.8))
    if pairs:
        data=np.array(pairs)
        axes[0].scatter(data[:,1],data[:,2],s=13,c=data[:,0],cmap='viridis',alpha=.5)
        lim=[max(1e-18,float(data[:,1:3].min())),float(data[:,1:3].max())];axes[0].plot(lim,lim,'k--',lw=1)
        axes[0].set_xscale('log');axes[0].set_yscale('log')
    axes[0].set_xlabel('Residual witness before enrichment');axes[0].set_ylabel('Same-direction witness after enrichment')
    if paths:
        data=np.array(paths);axes[1].scatter(data[:,0],data[:,1],s=15,alpha=.5)
        low=float(data[:,:2].min());high=float(data[:,:2].max());axes[1].plot([low,high],[low,high],'k--',lw=1)
    axes[1].set_xlabel('First-order finite-path prediction');axes[1].set_ylabel('Evaluated finite-path nonlinear bias B')
    for ax in axes:ax.grid(alpha=.2)
    fig.suptitle('Consequential directions: actual witnesses, not whole-class upper bounds');fig.tight_layout()
    for ext in ['png','svg']:fig.savefig(out/f'figure3_witness_nonlinearity.{ext}',dpi=160)
    plt.close(fig)
    dump(out/'figure_data.json',{'residual_pairs':pairs,'finite_path_rows':paths,'raw_results_sha256':hashlib.sha256((out/'raw_results.csv').read_bytes()).hexdigest()})
