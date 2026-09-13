from core import *
import re,sys,shutil
sys.path.insert(0,str(ROOT/'sirna_gnn_empirical/v1'))
from common import pairing,COMP
rows,raw=load();manifest=json.loads((OLD/'graph_manifest.json').read_text());vocab=manifest['modification_vocabulary'];source=OLD/'raw/davis2025_supplements/Supplemental Tables.xlsx';df=pd.read_excel(source,sheet_name='Supplemental Table 7',header=4);audit=pd.read_csv(RUN/'source_audit/davis_orientation_rows.csv');audit=audit[audit.sheet=='Supplemental Table 7'];assert len(audit)==118 and audit.guide_target_reverse_complement_13mer.all() and audit.guide_passenger_reverse_complement_13mer.all()
def km(s):return {s.replace('T','U')[i:i+13] for i in range(len(s)-12)}
known=set()
for x in rows:
 for st in ['guide','passenger']:known|=km(x[st]['sequence'])
new=[];checks=[];N=len(df);X=np.zeros((N,64,93),np.float32);A=np.zeros((N,8,64,64),np.float32);mask=np.zeros((N,64),np.float32)
for j,x in df.iterrows():
 rr=dict(record_id=f'davis:S7:{j+6}',dataset='DavisS7',source_family='Davis2025',split='Davis_S7',study_group='Davis2025_S7',activity=float(1-x.iloc[7]/100),measurement_sd=float(x.iloc[8]/100),dose_reported=0.,dose_unit='not_verified_for_this_row',time_h=None,cell=str(x.iloc[9]),target=str(x.iloc[1]),delivery='cholesterol-conjugated; passive uptake described in study',assay_id=None,raw_target=str(x.iloc[4]),accession=str(x.iloc[6]),compound=str(x.iloc[0]),raw_guide=str(x.iloc[2]),raw_passenger=str(x.iloc[3]))
 for st,(key,col) in enumerate([('guide',2),('passenger',3)]):
  literal=str(x.iloc[col]);tokens=list(re.finditer(r'\(([mfr])([ACGU])\)',literal));seq=''.join(t[2] for t in tokens);nodes=[];links=[];assert 14<=len(seq)<=32
  for i,t in enumerate(tokens):
   m={'m':'2-O-Methyl','f':'2-Fluoro','r':'unmodified_as_annotated'}[t[1]];nodes.append(dict(base=t[2],mods=[m],stereo='not_reported'));k=st*32+i;mask[j,k]=1;X[j,k,'ACGUT'.index(t[2])]=1;X[j,k,5+st]=1;X[j,k,7+i]=1;X[j,k,39+vocab.index(m)]=1;X[j,k,84]=1;X[j,k,87]=i==0;X[j,k,88]=i==len(tokens)-1
   if i<len(tokens)-1:
    between=literal[t.end():tokens[i+1].start()];assert between in ['', '#'];link='PS' if between=='#' else 'PO';links.append(link);rel=int(link=='PS');A[j,rel,k+1,k]=1;A[j,rel+3,k,k+1]=1
  if literal.startswith('P'):X[j,st*32,91]=1
  if literal.endswith('-TegChol'):X[j,st*32+len(tokens)-1,92]=1
  rr[key]=dict(sequence=seq,nodes=nodes,linkages=links,terminal=[]) # Terminal identities are retained in raw strings and graph flags.
 pp,fraction=pairing(rr['guide']['sequence'],rr['passenger']['sequence']);assert fraction>=.8
 for a,b in pp:
  rel=6 if COMP[rr['guide']['sequence'][a]]==rr['passenger']['sequence'][b] else 7;A[j,rel,a,32+b]=A[j,rel,32+b,a]=1
 rr['pairing_fraction']=fraction;rr['pairing']=pp;overlap=any(km(rr[st]['sequence'])&known for st in ['guide','passenger']);checks.append(dict(record_id=rr['record_id'],known_source_or_sequence_overlap=overlap,pairing_fraction=fraction));new.append(rr)
# Sequence grouping within the new cohort is independent of outcomes.
parent=list(range(N))
def root(i):
 while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
 return i
for i in range(N):
 for j in range(i):
  if any(km(new[i][a]['sequence'])&km(new[j][b]['sequence']) for a in ['guide','passenger'] for b in ['guide','passenger']):parent[root(i)]=root(j)
for i,x in enumerate(new):x['sequence_group']='davisS7_'+str(root(i))
keep=np.array([i for i,c in enumerate(checks) if not c['known_source_or_sequence_overlap']],int)
write(RUN/'external_admission.json',dict(source='Davis2025 Supplemental Table 7',source_workbook_sha256=sha(source),candidates=118,admitted=len(keep),known_13mer_overlap_excluded=118-len(keep),sequence_components=len({new[i]['sequence_group'] for i in keep}),no_pretrained_checkpoint=True,guide_conversion='None: literal printed order, primary m/f/r and # notation. All 118 guide strings and duplexes pass reverse-complement 13mer checks.',dose_time='Not assigned to S7 rows because main-text dose/time description does not unambiguously cover every S7 target. Both are masked by training-only context support; no within-assay ranking is reported.',endpoint='1 - reported native mRNA expression percent/100; no clipping; table SD/100 preserved, n unresolved for individual S7 rows.',status='Additional separately sourced observational endpoint evaluation with frozen all-ENsiRNA checkpoints. Public outcomes were available in an acquired workbook; this is not a prospective blind trial. No APP or S7 labels select models.',S1_status='quarantined',S6_status='24 rows fail orientation; 26 pass. S6 reporter outcomes are not mixed into the prespecified S7 native-outcome cohort.',B2='Not used for training: S7 is one named scaffold family and exact comparable dose/time/control dependency per row is unresolved. Original exploratory APP supervised protocol remains fixed.'))
pd.DataFrame(checks).to_csv(RUN/'external_overlap_checks.csv',index=False);(RUN/'external_observations.jsonl').write_text(''.join(json.dumps(x)+'\n' for x in new));np.savez_compressed(RUN/'external_graphs.npz',X=X,A=A,mask=mask)
allrows=rows+new;combined={k:np.concatenate([raw[k],{'X':X,'A':A,'mask':mask,'C':np.zeros((N,8),np.float32)}[k]]) for k in raw};tr=np.array([i for i,x in enumerate(rows) if x['dataset']=='ENsiRNA']);te=keep+len(rows);data,sup=transform(combined,allrows,tr,True,True);fits=json.loads((RUN/'deployment_fits.json').read_text());write(RUN/'external_freeze.json',dict(admission_sha256=sha(RUN/'external_admission.json'),input_sha256=sha(RUN/'external_graphs.npz'),checkpoints={f['checkpoint']:f['checkpoint_sha256'] for f in fits},metrics=['pooled MSE/MAE','equal sequence component MSE','per target/cell descriptive MSE','conditional sequence-component paired bootstrap'],selected_before_prediction=True));device=setup();predictions=[]
for f in fits:
 if f['kind']=='chemistry_tree':
  with (RUN/f['checkpoint']).open('rb') as ff:st=pickle.load(ff)
  degree=data['A'].sum(-1).copy();od=degree.copy();degree[:,2]=od[:,:3].sum(1);degree[:,5]=od[:,3:6].sum(1)
  for rel in [0,1,3,4,6,7]:degree[:,rel]*=sup['relation_support'][rel]
  xx=np.concatenate([data['X'].reshape(len(allrows),-1),degree.reshape(len(allrows),-1),data['C']],1);pp=st['model'].predict(xx[te][:,st['columns']])
 else:
  st=torch.load(RUN/f['checkpoint'],map_location=device,weights_only=False);m=model(f['kind'],True,st['support']['relation_support']).to(device);m.load_state_dict(st['model']);pp=pred(m,data,te,device,st['target_mean'],st['target_std'])
 for i,p in zip(te,pp):predictions.append(dict(record_id=allrows[i]['record_id'],model=f['kind'],seed=str(f['seed']),prediction=float(p)))
df=pd.DataFrame(predictions);df.to_csv(RUN/'external_predictions.csv',index=False);write(RUN/'external_prediction_commit.json',dict(predictions_sha256=sha(RUN/'external_predictions.csv')));en=df.groupby(['record_id','model'],as_index=False).prediction.mean();en['seed']='ensemble';df=pd.concat([df,en]);md=pd.DataFrame([{k:x[k] for k in ['record_id','activity','sequence_group','target','cell']} for x in new]);df=df.merge(md,on='record_id');df.to_csv(RUN/'external_predictions_with_labels.csv',index=False);metrics=[];details=[]
for (kind,seed),g in df.groupby(['model','seed']):
 e=(g.prediction-g.activity)**2;metrics.append(dict(model=kind,seed=seed,n=len(g),mse=e.mean(),mae=abs(g.prediction-g.activity).mean(),equal_component_mse=np.average(e,weights=weights(g.sequence_group))))
 for (t,c),a in g.groupby(['target','cell']):details.append(dict(model=kind,seed=seed,target=t,cell=c,n=len(a),mse=np.mean((a.prediction-a.activity)**2)))
pd.DataFrame(metrics).to_csv(RUN/'external_comparison.csv',index=False);pd.DataFrame(details).to_csv(RUN/'external_target_results.csv',index=False)
print(json.loads((RUN/'external_admission.json').read_text()));print(pd.DataFrame(metrics).query("seed=='ensemble'").to_string(index=False))
