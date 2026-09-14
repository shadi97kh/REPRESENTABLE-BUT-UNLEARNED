from core import *
df=pd.read_csv(RUN/'external_predictions_with_labels.csv');df=df[df.seed=='ensemble'];p=df.pivot(index='record_id',columns='model',values='prediction');md=df.drop_duplicates('record_id').set_index('record_id').loc[p.index];groups=sorted(md.sequence_group.unique());rng=np.random.default_rng(30260915);draws=rng.integers(0,len(groups),size=(2000,len(groups)));result=[]
for ref in ['chemistry_tree','token_cnn']:
 for kind in p.columns:
  if kind==ref:continue
  diff=(p[kind]-md.activity)**2-(p[ref]-md.activity)**2;nums=np.array([diff[md.sequence_group==g].sum() for g in groups]);den=np.array([(md.sequence_group==g).sum() for g in groups]);boot=nums[draws].sum(1)/den[draws].sum(1);lo,hi=np.quantile(boot,[.025,.975]);result.append(dict(model=kind,reference=ref,difference_mse=diff.mean(),lower95=lo,upper95=hi,observations=len(p),sequence_components=len(groups),bootstrap_replicates=2000))
pd.DataFrame(result).to_csv(RUN/'external_paired_intervals.csv',index=False);print(pd.DataFrame(result).to_string(index=False))
