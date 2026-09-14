"""Full support, fold, assay and primary-evidence exports from completed fits."""
from common import *
from engine import load
from metrics import metrics
rows,raw=load();meta=pd.DataFrame(rows);fitmap={};partitions={};fitparts=[];groups=[];feature=[];relations=[];param=[];roles=[]
for fp in sorted((RUN/'fits').rglob('fit.json')):
 r=json.loads(fp.read_text());fitmap[r['name']]=r;sig=r['signature'];tr=sig['train'];key=hashlib.sha256(json.dumps(tr).encode()).hexdigest()[:12];fitparts.append(dict(fit=r['name'],partition=key))
 if key not in partitions:
  partitions[key]=tr;sup=json.loads((fp.parent/'preprocessing.json').read_text());active=raw['X'][tr][raw['mask'][tr].astype(bool)]
  for j in range(93):feature.append(dict(partition=key,feature=j,active_training_nodes=int(len(active)),nonzero_training_nodes=int(np.count_nonzero(active[:,j])),retained=bool(sup['node_support'][j])))
  for g,z in meta.iloc[tr].groupby('study_group'):groups.append(dict(partition=key,study_group=g,rows=len(z),sequence_components=z.sequence_group.nunique(),source_labels='|'.join(sorted(set(z.source_family)))))
 if '/final/' in r['name']:
  for role in ['train','validation','test']:
   ii=sig[role];counts=raw['A'][ii].sum((0,2,3)) if ii else np.zeros(8);m=meta.iloc[ii]
   for j,c in enumerate(counts):relations.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],role=role,partition=key,relation=j,directed_edges=int(c),observations=len(ii),active_nodes=int(raw['mask'][ii].sum()) if ii else 0))
   roles.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],role=role,partition=key,observations=len(ii),study_linked_groups=m.study_group.nunique(),sequence_components=m.sequence_group.nunique(),source_families=m.source_family.nunique()))
  gp=fp.parent/'gradient_audit.json'
  if gp.exists() and r['method']!='token_cnn':
   gs=pd.DataFrame(json.loads(gp.read_text()));sup=json.loads((fp.parent/'preprocessing.json').read_text())
   # Relations with structurally unsupported transforms. Corrected base remains active when any directional bond exists.
   supported=np.array(sup['relation_support'],bool)
   for layer,z in gs.groupby('parameter'):
    unavailable=[i for i in range(8) if not supported[i]]
    if 'corrected' in r['method']:
     unavailable=[i for i in unavailable if i not in [2,5]]
     if not any(supported[:3]):unavailable.append(2)
     if not any(supported[3:6]):unavailable.append(5)
    param.append(dict(fit=r['name'],method=r['method'],seed=r['seed'],layer=layer,parameters_total=r['parameters'],unsupported_relation_scalars=1024*len(unavailable),zero_data_gradient_scalars=int(z.loc[z.data_gradient_norm==0,'scalars'].sum()),minimum_relation_gradient=float(z.data_gradient_norm.min()),maximum_relation_gradient=float(z.data_gradient_norm.max()),**{f'g{int(t.relation)}':float(t.data_gradient_norm) for t in z.itertuples()},**{f'delta{int(t.relation)}':float(t.parameter_l1_change) for t in z.itertuples()}))
export(pd.DataFrame(fitparts),'tables/fit_partition_mapping.csv');export(pd.DataFrame(groups),'tables/training_partition_groups.csv');export(pd.DataFrame(feature),'tables/feature_support_counts.csv');export(pd.DataFrame(relations),'tables/final_relation_counts.csv');export(pd.DataFrame(roles),'tables/final_role_counts.csv');export(pd.DataFrame(param),'tables/final_parameter_support.csv')
# Historical reference count explicitly separate from the revision's grouped training populations.
hist=np.load(OLD/'train_graphs.npz');hc=hist['A'].sum((0,2,3));export(pd.DataFrame(dict(relation=range(8),directed_edges=hc.astype(int))),'tables/historical_training_relation_counts.csv')
ens=pd.read_csv(RUN/'evaluated/activity_ensemble_predictions.csv');seeds=pd.read_csv(RUN/'evaluated/activity_seed_predictions.csv');assays=[];folds=[]
for (cohort,method,assay),z in ens.groupby(['cohort','method','assay_id'],dropna=False):
 assays.append(dict(cohort=cohort,method=method,assay=str(assay),sequence_components=z.sequence_group.nunique(),study_linked_groups=z.study_group.nunique(),source_families=z.source_family.nunique(),**metrics(z.activity,z.prediction)))
# Per-assay oracle references are descriptive; fitted constants retain their actual training means.
for (cohort,assay),z in ens.drop_duplicates(['cohort','record_id']).groupby(['cohort','assay_id'],dropna=False):assays.append(dict(cohort=cohort,method='oracle_test_mean_diagnostic',assay=str(assay),sequence_components=z.sequence_group.nunique(),study_linked_groups=z.study_group.nunique(),source_families=z.source_family.nunique(),**metrics(z.activity,np.repeat(z.activity.mean(),len(z)))))
for (cohort,method,pop,seed),z in seeds.groupby(['cohort','method','population','seed']):folds.append(dict(cohort=cohort,method=method,population=pop,seed=seed,**metrics(z.activity,z.prediction)))
for (cohort,pop),z in seeds.drop_duplicates(['cohort','population','record_id']).groupby(['cohort','population']):
 for field in ['training_row_mean','training_equal_group_mean']:
  folds.append(dict(cohort=cohort,method=field,population=pop,seed='training_constant',**metrics(z.activity,z[field])))
 folds.append(dict(cohort=cohort,method='oracle_test_mean_diagnostic',population=pop,seed='oracle',**metrics(z.activity,np.repeat(z.activity.mean(),len(z)))))
export(pd.DataFrame(assays),'tables/per_assay_metrics.csv');export(pd.DataFrame(folds),'tables/per_fold_seed_metrics.csv')
# Selection records: every configuration score, selected configurations and pair-loss weights.
selection=[]
for stage in ['activity','pair']:
 for q in json.loads((RUN/f'{stage}_selection.json').read_text()):
  for cfg,score in q.get('scores',q.get('inner_scores',{})).items():selection.append(dict(stage=stage,population=q['population'],method=q['method'],kind=q.get('kind','activity_config'),configuration=cfg,validation_score=score,selected=str(cfg)==str(q.get('config',{}).get('name',q.get('lambda_pair')))))
export(pd.DataFrame(selection),'tables/selection_scores.csv')
# Audit all fixed final-seed assignments, optimizer records, and output identities.
f=pd.read_csv(RUN/'tables/all_fit_summary.csv');neural=f[f.executed_updates>0];write(RUN/'tables/optimization_checks.json',dict(fits=len(f),neural_fits=len(neural),neural_fits_without_updates=int(((f.method.isin(['original_gnn','corrected_gnn','original_no_message','corrected_no_message','token_cnn']))&(f.executed_updates<=0)).sum()),selected_le_executed=bool((neural.selected_updates<=neural.executed_updates).all()),selected_epoch_le_executed=bool((neural.selected_epoch<=neural.epochs_executed).all()),minimum_executed_epochs=int(neural.epochs_executed.min()),maximum_executed_epochs=int(neural.epochs_executed.max()),finite_unfavorable_seeds_discarded=0,final_seed_schedule=SEEDS))
# Primary panel states and complete observed rectangle values are typed table inputs.
states=json.loads((RUN/'b3/primary_strand_states.json').read_text());used=readlines(RUN/'b3_primary/observations.jsonl');names=sorted({r[k] for r in used for k in ['AS','SS']});st=[]
for name in names:
 q=states[name];st.append(dict(state=name,strand=q['strand'],sequence=q['sequence'],modifications='; '.join(f"{i+1}:"+'+'.join(n['mods']) for i,n in enumerate(q['nodes']) if n['mods']) or 'unmodified'))
export(pd.DataFrame(st),'tables/b3_primary_states.csv');rect=json.loads((RUN/'b3_primary/rectangles.json').read_text());export(pd.DataFrame([dict(rectangle_id=q['rectangle_id'],AS=q['AS'],SS=q['SS'],**{f'y{k}':v for k,v in q['observed'].items()},interaction=q['observed_interaction'],**{f'sd{k}':v for k,v in q['endpoint_sd'].items()}) for q in rect]),'tables/b3_primary_measurements.csv')
write(RUN/'analysis_supplement_summary.json',dict(unique_training_partitions=len(partitions),fit_partition_links=len(fitparts),final_relation_rows=len(relations),feature_count_rows=len(feature),per_assay_rows=len(assays),per_fold_seed_rows=len(folds),all_main_numbers_must_be_generated=True))
print('Supplement complete',len(partitions),len(assays),len(folds))

structural=[]
for name,r in fitmap.items():
 if '/final/' not in name or r['method'] not in ['original_gnn','corrected_gnn','original_no_message','corrected_no_message']:continue
 pp=RUN/'fits'/name;sup=json.loads((pp/'preprocessing.json').read_text());tr=r['signature']['train'];unused_slots=int((raw['mask'][tr].sum(0)==0).sum());inactive_input=32*sum(not x for x in sup['node_support']);inactive_readout=64*(32*unused_slots+sum(not x for x in sup['context_varying']));unsupported=sum(z['unsupported_relation_scalars'] for z in param if z['fit']==name);inactive=inactive_input+inactive_readout+unsupported
 structural.append(dict(fit=name,method=r['method'],seed=r['seed'],total_scalars=r['parameters'],inactive_input_scalars=inactive_input,inactive_readout_scalars=inactive_readout,inactive_relation_scalars=unsupported,structurally_inactive_scalars=inactive,train_reachable_upper_bound=r['parameters']-inactive,qualification='Upper bound on data-reachable scalars, not proof every remaining scalar has nonzero realized gradient'))
export(pd.DataFrame(structural),'tables/structural_parameter_counts.csv')
bp=pd.read_csv(RUN/'pair_predictions.csv');seedcomp=[]
for a,b in [('pair_gnn','training_mean'),('pair_gnn','pair_no_message'),('pair_gnn','pair_tree'),('pair_gnn','activity_gnn')]:
 for seed,z in bp[bp.method==a].groupby('seed'):
  zz=bp[(bp.method==b)&(bp.seed==(0 if b=='training_mean' else seed))];q=z.merge(zz[['pair_id','prediction']],on='pair_id',suffixes=('_a','_b'),validate='one_to_one');d=(q.prediction_a-q.observed_difference)**2-(q.prediction_b-q.observed_difference)**2
  for label,w in [('rows',np.ones(len(q))),('equal_sequence',weights(q.sequence_group))]:seedcomp.append(dict(a=a,b=b,seed=seed,weighting=label,difference=float(np.average(d,weights=w))))
export(pd.DataFrame(seedcomp),'tables/pair_seed_differences.csv')
# Compact original negative-result reproduction, without mutating historical files.
hp=pd.read_csv(OLD/'heldout_predictions.csv').merge(meta[['record_id','activity','sequence_group']],on='record_id',validate='many_to_one');history=[]
for (split,model),z in hp.groupby(['split','model']):history.append(dict(version='v1_frozen',cohort=split,method=model,estimand='activity',weighting='rows',**metrics(z.activity,z.prediction)))
hb=pd.read_csv(OLD/'chemical_effect_metrics.csv')
for z in hb.to_dict('records'):history.append(dict(version='v1_frozen',cohort='APP_B2',method=z['model'],estimand='measured_bundle_difference',weighting='rows',mse=z['pair_mse'],n=z['pairs']))
export(pd.DataFrame(history),'tables/historical_v1_results.csv')
write(RUN/'historical_v1_recomputation.json',dict(prediction_source=str((OLD/'heldout_predictions.csv').relative_to(ROOT)),prediction_sha256=sha(OLD/'heldout_predictions.csv'),B2_metric_source=str((OLD/'chemical_effect_metrics.csv').relative_to(ROOT)),B2_metric_sha256=sha(OLD/'chemical_effect_metrics.csv'),interpretation='Activity metrics recomputed from preserved ensemble predictions. B2 errors copied from exact preserved metric table; not a new fit or revised B2 protocol.'))
# Within-assay variation diagnostic: remove shared assay means only for this labeled diagnostic.
# It is not deployed recalibration and never enters training/model selection.
predpair=pd.read_csv(RUN/'evaluated/pair_ensemble_predictions.csv').merge(meta[['record_id','assay_id']],left_on='reference_id',right_on='record_id',validate='many_to_one');within=[]
for method,z in predpair.groupby('method'):
 if method=='training_majority':continue
 y=z.observed_difference.to_numpy(float);p=z.prediction.to_numpy(float);_,aa=np.unique(z.assay_id.astype(str),return_inverse=True);unique,gg=np.unique(z.sequence_group,return_inverse=True);base=np.ones(len(z));rng=np.random.default_rng(90413)
 def centered(w):
  total=np.bincount(aa,weights=w);ym=np.divide(np.bincount(aa,weights=w*y),total,out=np.zeros(len(total)),where=total>0);pm=np.divide(np.bincount(aa,weights=w*p),total,out=np.zeros(len(total)),where=total>0);return metrics(y-ym[aa],p-pm[aa],w)
 point=centered(base);corr=[]
 for _ in range(10000):
  multiplicity=np.bincount(rng.integers(len(unique),size=len(unique)),minlength=len(unique));v=centered(multiplicity[gg])['correlation']
  if v is not None:corr.append(v)
 within.append(dict(method=method,assays=len(set(aa)),sequence_components=len(unique),**point,conditional_correlation_lower=float(np.quantile(corr,.025)) if corr else None,conditional_correlation_upper=float(np.quantile(corr,.975)) if corr else None,finite_bootstrap_correlations=len(corr),diagnostic='Both measured and predicted assay means removed within each resample; no fitted or deployable intercept correction',analysis_seed=90413,bootstrap_resamples=10000))
export(pd.DataFrame(within),'tables/B2_within_assay_variation.csv')
# The controlled support rule can change held-out routing while leaving ENsiRNA fitting identical.
# Inspect saved states/predictions; no refit or new inference is needed.
matched=[]
for name,r in fitmap.items():
 if not name.startswith('activity/') or '/final/' not in name or r['method'] not in ['original_gnn','original_no_message']:continue
 other=name.replace('original_gnn','corrected_gnn').replace('original_no_message','corrected_no_message');q=fitmap[other];p0=RUN/'fits'/name;p1=RUN/'fits'/other;v0=np.load(p0/'predictions.npz')['validation'];v1=np.load(p1/'predictions.npz')['validation'];samecfg=r['config']==q['config'];delta=None;identical=None
 if r['status']==q['status']=='completed':
  import torch
  a=torch.load(p0/'best.pt',map_location='cpu',weights_only=False)['model'];b=torch.load(p1/'best.pt',map_location='cpu',weights_only=False)['model'];delta=max(float((a[k]-b[k]).abs().max()) for k in a);identical=all(torch.equal(a[k],b[k]) for k in a)
 matched.append(dict(original_fit=name,corrected_fit=other,population=name.split('/')[1],seed=r['seed'],same_training_and_validation=r['signature']['train']==q['signature']['train'] and r['signature']['validation']==q['signature']['validation'],same_selected_configuration=samecfg,same_selected_epoch=r['selected_epoch']==q['selected_epoch'],same_executed_updates=r['executed_updates']==q['executed_updates'],validation_prediction_max_absolute_difference=float(np.max(abs(v0-v1))),selected_parameter_max_absolute_difference=delta,selected_parameters_bit_identical=identical,qualification='Saved-state comparison only. Different inference routing can change APP/S7 predictions even with identical selected tensors.'))
export(pd.DataFrame(matched),'tables/matched_support_training_identity.csv')

# Complete held-out pair-fold diagnostics, keeping pooled OOF statistics separate.
bfold=[];efold=[]
def pairfold(z,method,seed,pop):
 observed=z.observed_difference.to_numpy(float);predicted=z.prediction.to_numpy(float);obs=np.where(abs(observed)<=.02,0,np.sign(observed));calls=np.where(abs(predicted)<=.02,0,np.sign(predicted))
 for weighting,key in [('rows',None),('equal_sequence','sequence_group'),('equal_background','background_id')]:
  w=np.ones(len(z)) if key is None else weights(z[key]);q=metrics(observed,predicted,w)
  if method=='training_majority':
   for field in ['mse','mae','r_squared','correlation','prediction_sd','mean_bias','prediction_mean']:q[field]=None
  bfold.append(dict(population=pop,method=method,seed=seed,weighting=weighting,sequence_components=z.sequence_group.nunique(),exact_backgrounds=z.background_id.nunique(),**q,sign_accuracy=float(np.average(obs==calls,weights=w)),non_tie_correct=int(((obs==calls)&(obs!=0)).sum()),non_ties=int((obs!=0).sum()),negative_calls=int((calls<0).sum()),zero_calls=int((calls==0).sum()),positive_calls=int((calls>0).sum())))
for (pop,method,seed),z in bp.groupby(['population','method','seed']):pairfold(z,method,seed,pop)
for (pop,method),z in pd.read_csv(RUN/'evaluated/pair_ensemble_predictions.csv').groupby(['population','method']):pairfold(z,method,'ensemble',pop)
ep=pd.read_csv(RUN/'pair_endpoint_predictions.csv')
for (pop,method,seed),z in ep.groupby(['population','method','seed']):
 efold.append(dict(population=pop,method=method,seed=seed,sequence_components=z.sequence_group.nunique(),**metrics(z.activity,z.prediction,weights(z.sequence_group))))
for (pop,method),z in ep.groupby(['population','method']):
 q=z.groupby('record_id',as_index=False).agg(prediction=('prediction','mean'),activity=('activity','first'),sequence_group=('sequence_group','first'));efold.append(dict(population=pop,method=method,seed='ensemble',sequence_components=q.sequence_group.nunique(),**metrics(q.activity,q.prediction,weights(q.sequence_group))))
export(pd.DataFrame(bfold),'tables/B2_fold_seed_metrics.csv');export(pd.DataFrame(efold),'tables/B2_endpoint_fold_metrics.csv')
# Closed B2 identity dictionaries make measured endpoints assessable in the PDF too.
lookup={r['record_id']:r for r in rows};original_pairs=readlines(OLD/'eligible_pairs.jsonl');templates={};duplexes={};endpoint_code={};template_records=[];duplex_records=[];pair_records=[]
for rid in sorted({p[k] for p in original_pairs for k in ['reference_id','changed_id']}):
 row=lookup[rid];states=[]
 for strand in ['guide','passenger']:
  st=row[strand];key=json.dumps(dict(nodes=[dict(mods=n['mods'],stereo=n['stereo']) for n in st['nodes']],linkages=st['linkages'],terminal=st['terminal']),sort_keys=True)
  if key not in templates:
   code='C'+str(len(templates)+1);templates[key]=code;template_records.append(dict(template=code,length=len(st['nodes']),positional_chemistry='; '.join(str(i+1)+': '+', '.join(n['mods'])+(' [stereo '+n['stereo']+']' if n['stereo']!='not_reported' else '') for i,n in enumerate(st['nodes'])),linkages='; '.join(str(i+1)+'-'+str(i+2)+': '+v for i,v in enumerate(st['linkages'])),terminals=json.dumps(st['terminal'],sort_keys=True),stereochemistry='Unreported except explicit brackets'))
  states.append(templates[key])
 key=(row['guide']['sequence'],row['passenger']['sequence'],*states)
 if key not in duplexes:
  code='D'+str(len(duplexes)+1);duplexes[key]=code;duplex_records.append(dict(duplex=code,guide_sequence=key[0],passenger_sequence=key[1],guide_chemistry=key[2],passenger_chemistry=key[3],sequence_group=row['sequence_group']))
 endpoint_code[rid]=duplexes[key]
for p in original_pairs:
 ref=lookup[p['reference_id']];changed=lookup[p['changed_id']];assert ref['guide']['sequence']==changed['guide']['sequence'] and ref['passenger']['sequence']==changed['passenger']['sequence'];pair_records.append(dict(pair_id=p['pair_id'],background_id=p['background_id'],reference_id=p['reference_id'],changed_id=p['changed_id'],reference_duplex=endpoint_code[p['reference_id']],changed_duplex=endpoint_code[p['changed_id']],assay_id=p['assay_id'],reference_activity=ref['activity'],changed_activity=changed['activity']))
export(pd.DataFrame(template_records),'tables/B2_chemistry_templates.csv');export(pd.DataFrame(duplex_records),'tables/B2_duplex_states.csv');export(pd.DataFrame(pair_records),'tables/B2_measured_endpoint_mapping.csv')
