"""Frozen post-preprocessing identities; collisions are not inferred from predictions."""
from engine import *
from metrics import metrics
import collections
V3=ROOT/'runs/sirna_gnn_empirical/v3-20260914T013314Z';out=RUN/'visibility';out.mkdir(exist_ok=True);assert (RUN/'fit_completion.json').exists()
rows,raw=load();primary=readlines(V3/'b3_primary/observations.jsonl');praw=dict(np.load(V3/'b3_primary/graphs.npz'));both=rows+primary;n=len(rows);data={k:np.concatenate([raw[k],praw[k]]) for k in raw};bidx=np.arange(n,len(both));sidx=np.array([i for i,x in enumerate(rows) if x['record_id'].startswith('davis:')]);rect=json.loads((V3/'b3_primary/rectangles.json').read_text());oldpred=pd.read_csv(V3/'b3_primary/interaction_predictions.csv');oldpred=oldpred[oldpred.seed.astype(str)=='ensemble'];methods=['corrected_gnn','corrected_no_message','chemistry_tree','token_cnn'];vocab=json.loads((OLD/'graph_manifest.json').read_text())['modification_vocabulary'];summary=[];identities=[];forced=[];chem=[];metrics_out=[];leave=[];reference=[]
canonical_seen={}
def digest_arrays(parts):
 h=hashlib.sha256();chunks=[]
 for name,value in sorted(parts.items()):
  a=np.asarray(value,dtype='<f4').copy();a[a==0]=0;header=json.dumps([name,list(a.shape),'little-endian float32; signed zero canonicalized'],separators=(',',':')).encode();payload=a.tobytes(order='C');h.update(header);h.update(payload);chunks.extend([header,payload])
 key=h.hexdigest();canonical=tuple(chunks)
 if key in canonical_seen:assert canonical_seen[key]==canonical,'SHA256 collision: exact identity not established'
 else:canonical_seen[key]=canonical
 return key
cases=[('B3_original_source_heldout',V3/'fits/activity/outer0/final',bidx),('S7_original_deployment',V3/'fits/activity/deployment/final',sidx),('S7_source_excluded_reused',V3/'fits/activity/outer0/final',sidx),('S7_primary_reanchored',RUN/'fits/primary_reanchored/deployment/final',sidx)]
for case,base,idx in cases:
 for method in methods:
  canonical_seen.clear()
  files=[base/method/f's{seed}/fit.json' for seed in SEEDS];recs=[json.loads(p.read_text()) for p in files];assert len({json.dumps(r['signature']['train']) for r in recs})==1
  tr=np.array(recs[0]['signature']['train']);dd,sup=transform(data,both,tr)
  for fp in files:assert sup==json.loads((fp.parent/'preprocessing.json').read_text())
  if method=='chemistry_tree':
   columns=[]
   for fp in files:
    with (fp.parent/'model.pkl').open('rb') as f:st=pickle.load(f)
    columns.append(st['columns'])
   assert all(np.array_equal(columns[0],c) for c in columns);features=classical_features(dd,method)[idx][:,columns[0]]
  mapping={};route_mapping={}
  for j,i in enumerate(idx):
   X=dd['X'][i];A=dd['A'][i];mask=dd['mask'][i];C=dd['C'][i];carrier=dict(X=X,mask=mask,C=C)
   if method=='chemistry_tree':carrier=dict(features=features[j])
   elif method=='corrected_no_message':carrier['relation_receiver_degrees']=A.sum(-1)
   elif method=='corrected_gnn':carrier['A']=A
   actual=digest_arrays(carrier);routed=dict(carrier)
   if method in ['corrected_gnn','corrected_no_message']:
    B=A.copy();support=np.array(sup['relation_support'],float);B[2]=A[0]+A[1]+A[2];B[5]=A[3]+A[4]+A[5]
    for r in [0,1,3,4,6,7]:B[r]*=support[r]
    routed.pop('A',None);routed.pop('relation_receiver_degrees',None);routed['degree_denominator']=np.maximum(A.sum((0,2)),1);routed['effective_relations']=B if method=='corrected_gnn' else B.sum(-1)
   route=digest_arrays(routed);x=both[i];mapping[x['record_id']]=actual;route_mapping[x['record_id']]=route;identities.append(dict(case=case,method=method,record_id=x['record_id'],input_sha256=actual,routed_real_arithmetic_sha256=route,activity=x['activity']))
  groups=collections.defaultdict(list)
  for i in idx:groups[mapping[both[i]['record_id']]].append(float(both[i]['activity']))
  floor=sum(sum((np.array(v)-np.mean(v))**2) for v in groups.values())/len(idx)
  summary.append(dict(case=case,method=method,rows=len(idx),exact_inputs=len(groups),routed_real_arithmetic_inputs=len(set(route_mapping.values())),collision_classes=sum(len(v)>1 for v in groups.values()),rows_in_collision_classes=sum(len(v) for v in groups.values() if len(v)>1),largest_collision_class=max(map(len,groups.values())),finite_sample_collision_MSE_floor=floor,checked_frozen_preprocessors=10))
  for m in vocab:
   nodes=sum(m in node['mods'] for i in idx for st in ['guide','passenger'] for node in both[i][st]['nodes'])
   if nodes:chem.append(dict(case=case,method=method,chemical_name=m,observed_nodes=nodes,named_feature_supported=bool(sup['node_support'][39+vocab.index(m)]),qualification='Unsupported type name can disappear while its position/base/unmodified-indicator difference remains visible.'))
  if case.startswith('B3'):
   pr=oldpred[oldpred.method==method].set_index('rectangle_id');local=[]
   for r in rect:
    h={k:mapping[v] for k,v in r['record_ids'].items()};rh={k:route_mapping[v] for k,v in r['record_ids'].items()};f=collections.Counter([h['11'],h['00']])==collections.Counter([h['10'],h['01']]);rf=collections.Counter([rh['11'],rh['00']])==collections.Counter([rh['10'],rh['01']]);z=dict(method=method,rectangle_id=r['rectangle_id'],AS=r['AS'],SS=r['SS'],distinct_exact_inputs=len(set(h.values())),AS_marginal_invisible=h['10']==h['00'],SS_marginal_invisible=h['01']==h['00'],forced_exact_input_cancellation=f,forced_routed_real_arithmetic_cancellation=rf,observed_interaction=r['observed_interaction'],predicted_interaction=float(pr.loc[r['rectangle_id'],'prediction']));local.append(z);forced.append(z)
   frame=pd.DataFrame(local)
   for label,q in [('all',frame),('exact_forced',frame[frame.forced_exact_input_cancellation]),('not_exact_forced',frame[~frame.forced_exact_input_cancellation])]:
    if len(q):
     for pred in ['model','zero']:
      metrics_out.append(dict(method=method,stratum=label,predictor=pred,**metrics(q.observed_interaction,q.predicted_interaction if pred=='model' else np.zeros(len(q)))))
   for axis in ['AS','SS']:
    for value in sorted(frame[axis].unique()):
     q=frame[frame[axis]!=value];leave.append(dict(method=method,omitted_axis=axis,omitted_state=value,n=len(q),model_mse=np.mean((q.predicted_interaction-q.observed_interaction)**2),zero_mse=np.mean(q.observed_interaction**2),difference=np.mean((q.predicted_interaction-q.observed_interaction)**2-q.observed_interaction**2),interpretation='dependent descriptive leave-one-margin sensitivity; not independent folds or replicates'))
   sd=float(rect[0]['endpoint_sd']['00']);p=frame.predicted_interaction.to_numpy();y=frame.observed_interaction.to_numpy()
   for multiplier in [-1,0,1]:
    shift=multiplier*sd;yy=y+shift;reference.append(dict(method=method,reference_shift=shift,multiplier_of_reported_reference_SD=multiplier,model_mse=np.mean((p-yy)**2),zero_mse=np.mean(yy**2),model_minus_zero=np.mean((p-yy)**2-yy**2),interpretation='descriptive common-reference perturbation, not a confidence interval; unknown mean SE/covariance'))
  del dd
for name,values in [('input_counts',summary),('state_hashes',identities),('chemical_feature_visibility',chem),('B3_rectangle_visibility',forced),('B3_stratified_metrics',metrics_out),('B3_leave_one_margin',leave),('B3_shared_reference_sensitivity',reference)]:export(pd.DataFrame(values),f'visibility/{name}.csv')
write(out/'conventions.json',dict(hash_representation='Named arrays, shapes, float32 little-endian C-order; signed zeros normalized; no approximate-prediction collision rule',actual_carrier='GNN: X,A,mask,C. No-message: X, exact receiver relation counts,mask,C. CNN: X,mask,C. Tree: exact retained feature columns. Frozen training-only feature/context support applied.',routed_hash='An additional real-arithmetic quotient replaces equal backbone transforms by shared bases and gates unsupported residuals; original total degree remains. Algebraic regrouping may differ by floating-point roundoff; exact-input and algebraic flags reported separately.',forced_cancellation='Equality of multisets {11,00} and {10,01} implies zero signed contrast for any deterministic function of those exact inputs. Distinct inputs do not prove useful expressivity or learned chemistry.',collision_floor='Within-class empirical response variance minimizes squared error over deterministic functions of the observed encoded inputs; finite-sample identity, no population lower-bound claim.',B3='Original fixed-checkpoint predictions; 140 distinct dependent contrasts, 14 AS replacements, 10 SS replacements, one background and common reference. No independent-rectangle bootstrap.',new_fits=0))
print(pd.DataFrame(summary).to_string(index=False));print(pd.DataFrame(forced).groupby('method')[['forced_exact_input_cancellation','forced_routed_real_arithmetic_cancellation']].sum().to_string())
