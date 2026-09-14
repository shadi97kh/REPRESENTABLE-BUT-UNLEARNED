"""Label-independent graph construction. Fixed chemical vocabulary, training-only scaling."""
import argparse,json,collections,math,itertools
from pathlib import Path
import numpy as np
from scipy import sparse
from common import *
ap=argparse.ArgumentParser();ap.add_argument('--run',type=Path,required=True);r=ap.parse_args().run
SPLITS=['train','validation','test_internal','test_APP'];data={s:read_jsonl(r/f'{s}_inputs.jsonl') for s in SPLITS}
L=32;N=2*L;BASE='ACGUT';names=sorted({m for rows in data.values() for x in rows for strand in ['guide','passenger'] for n in x[strand]['nodes'] for m in n['mods']});vocab={m:i for i,m in enumerate(names)}
# Public, label-free input-schema dictionary preserves every named modification. It is not fitted to efficacy.
base_end=5;strand_end=7;pos_end=39;mod_end=pos_end+len(names);F=mod_end+9
context_names=['log1p_reported_concentration','concentration_units_unverified','time_h_div24','time_missing','context_metadata_available','pairing_fraction','guide_length_div32','passenger_length_div32']
contexts={}
for split,rows in data.items():
 contexts[split]=np.array([[math.log1p(x['dose_reported']),float(x['dose_unit']!='nM'),(x['time_h'] or 0)/24,float(x['time_h'] is None),float(x['cell'] is not None),x['pairing_fraction'],len(x['guide']['sequence'])/32,len(x['passenger']['sequence'])/32] for x in rows],dtype=np.float32)
mu=contexts['train'].mean(0);sd=contexts['train'].std(0);sd[sd<1e-6]=1
outputs=[];stats={};train_seen=set()
for x in data['train']:
 for st in ['guide','passenger']:
  for n in x[st]['nodes']:train_seen.update(n['mods'])
for split,rows in data.items():
 X=np.zeros((len(rows),N,F),np.float32);A=np.zeros((len(rows),8,N,N),np.float32);mask=np.zeros((len(rows),N),np.float32)
 guide_sparse=[];pair_sparse=[];unseen_rows=0
 for b,x in enumerate(rows):
  unseen=False
  for st,strand in enumerate(['guide','passenger']):
   state=x[strand];assert len(state['nodes'])<=L
   for i,node in enumerate(state['nodes']):
    j=st*L+i;mask[b,j]=1;X[b,j,BASE.index(node['base'])]=1;X[b,j,5+st]=1;X[b,j,7+i]=1
    for m in node['mods']:X[b,j,pos_end+vocab[m]]=1;unseen|=m not in train_seen
    X[b,j,mod_end+0]=node['stereo']=='not_reported';X[b,j,mod_end+1]=node['stereo']=='S_GNA';X[b,j,mod_end+2]=bool(node.get('linkage_direction_unresolved'));X[b,j,mod_end+3]=i==0;X[b,j,mod_end+4]=i==len(state['nodes'])-1
   for i,link in enumerate(state['linkages']):
    rel={'PO':0,'PS':1,'PO_assumed_from_annotation':2}[link];a=st*L+i;c=a+1;A[b,rel,c,a]=1;A[b,3+rel,a,c]=1
   for i,term in state['terminal']:
    j=st*L+min(i,len(state['nodes'])-1)
    if 'L96' in term:X[b,j,mod_end+5]=1
    elif 'vinyl' in term:X[b,j,mod_end+6]=1
    elif term=='5-Phosphate':X[b,j,mod_end+7]=1
    else:X[b,j,mod_end+8]=1
  for i,j in x['pairing']:
   a=i;c=L+j;rel=6 if COMP[x['guide']['sequence'][i]]==x['passenger']['sequence'][j].replace('T','U') else 7;A[b,rel,a,c]=A[b,rel,c,a]=1
  gs=x['guide']['sequence'];features=[i*5+BASE.index(v) for i,v in enumerate(gs)]
  guide_sparse.append(features)
  ps=features.copy()
  for i,j in itertools.combinations(range(min(24,len(gs))),2):
   pairindex=i*(47-i)//2+(j-i-1);ps.append(L*5+pairindex*25+BASE.index(gs[i])*5+BASE.index(gs[j]))
  pair_sparse.append(ps);unseen_rows+=unseen
 C=(contexts[split]-mu)/sd
 # Reference ridge receives guide bases plus the same observed numeric context; explicit chemistry omission.
 for name,inds,dim in [('guide',guide_sparse,L*5),('pairwise',pair_sparse,L*5+276*25)]:
  rr=[];cc=[]
  for i,z in enumerate(inds):rr.extend([i]*len(z));cc.extend(z)
  mat=sparse.csr_matrix((np.ones(len(rr),np.float32),(rr,cc)),shape=(len(rows),dim));mat=sparse.hstack([mat,sparse.csr_matrix(C)],format='csr');p=r/f'{split}_{name}.npz';sparse.save_npz(p,mat);outputs.append(p.name)
 p=r/f'{split}_graphs.npz';np.savez_compressed(p,X=X,A=A,mask=mask,C=C);outputs.append(p.name)
 write_json(r/f'{split}_graph_ids.json',[x['record_id'] for x in rows]);outputs.append(f'{split}_graph_ids.json')
 stats[split]={'observations':len(rows),'nodes':int(mask.sum()),'directed_edges':int(A.sum()),'rows_with_unseen_modification_names':int(unseen_rows),'bytes_uncompressed':X.nbytes+A.nbytes+mask.nbytes+C.nbytes}
 if split in ['train','validation']:
  ys=read_jsonl(r/f'{split}_observations.jsonl');assert [x['record_id'] for x in ys]==[x['record_id'] for x in rows]
  cnt=collections.Counter(x['study_group'] for x in ys);w=np.array([1/cnt[x['study_group']] for x in ys],np.float32);w*=len(w)/w.sum()
  p=r/f'{split}_targets.npz';np.savez(p,y=np.array([x['activity'] for x in ys],np.float32),weight=w);outputs.append(p.name)
manifest={'max_nodes':N,'max_strand_length':L,'node_features':F,'base_features':list(BASE),'position_columns':[7,pos_end],'chemistry_columns':[pos_end,mod_end],'modification_vocabulary':names,'vocabulary_exposure':'Fixed label-free schema includes all admitted input modification names. Test efficacy never used. Unseen names remain distinct representation columns; effective training support is reported.','train_seen_modifications':sorted(train_seen),'context_features':context_names,'context_mean_training':mu.tolist(),'context_sd_training':sd.tolist(),'edge_relations':['backbone_forward_PO','backbone_forward_PS','backbone_forward_unresolved','backbone_reverse_PO','backbone_reverse_PS','backbone_reverse_unresolved','antiparallel_canonical','antiparallel_mismatch'],'pairing_algorithm':'Choose overlap>=14 maximizing number of Watson-Crick matches, then minimizing mismatches, maximizing overlap, minimizing absolute offset and offset. Admit match fraction>=0.8. DNA T pairs with A, remains distinct from U in features. Edges are inferred, not measured.','terminal_columns':[mod_end+5,F],'metadata_columns':[mod_end,mod_end+5],'statistics':stats,'test_targets_in_graph_files':False}
write_json(r/'graph_manifest.json',manifest);outputs.append('graph_manifest.json');write_json(r/'graphs.outputs.json',outputs);print(json.dumps(manifest,indent=2))
