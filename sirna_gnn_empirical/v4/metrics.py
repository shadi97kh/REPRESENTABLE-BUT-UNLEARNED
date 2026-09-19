"""Declared estimands; fixed-prediction component bootstrap and exact fractional ties."""
from common import *
from scipy.stats import rankdata
ANALYSIS_SEED=90413
BOOTSTRAPS=10000
def metrics(y,p,w=None):
 y=np.asarray(y,float);p=np.asarray(p,float);w=np.ones(len(y)) if w is None else np.asarray(w,float);w=w/w.sum()
 ym=w@y;pm=w@p;vy=w@((y-ym)**2);vp=w@((p-pm)**2);mse=w@((p-y)**2)
 return dict(n=len(y),mse=mse,mae=w@abs(p-y),r_squared=1-mse/vy if vy>1e-20 else None,correlation=(w@((y-ym)*(p-pm)))/np.sqrt(vy*vp) if vy*vp>1e-20 else None,prediction_sd=np.sqrt(vp),label_sd=np.sqrt(vy),mean_bias=pm-ym,label_mean=ym,prediction_mean=pm,oracle_test_mean_mse=vy)
def bootstrap_difference(df,a,b,group='sequence_group',weighting='rows'):
 z=df[['record_id',group,'activity',a,b]].copy();z['loss_difference']=(z[a]-z.activity)**2-(z[b]-z.activity)**2
 sums=z.groupby(group).loss_difference.agg(['sum','count','mean']);rng=np.random.default_rng(ANALYSIS_SEED);n=len(sums);val=[]
 # Sampling multiplicities preserve all records of each sequence component.
 for start in range(0,BOOTSTRAPS,100):
  ix=rng.integers(n,size=(min(100,BOOTSTRAPS-start),n))
  v=sums['sum'].to_numpy()[ix].sum(1)/sums['count'].to_numpy()[ix].sum(1) if weighting=='rows' else sums['mean'].to_numpy()[ix].mean(1);val.extend(v)
 est=z.loss_difference.mean() if weighting=='rows' else sums['mean'].mean()
 return dict(a=a,b=b,difference=est,lower=np.quantile(val,.025),upper=np.quantile(val,.975),groups=n,weighting=weighting,bootstrap_resamples=BOOTSTRAPS,analysis_seed=ANALYSIS_SEED,uncertainty='conditional on fitted predictions and observed cohort'),z

def ranking(df,pred='prediction',k=5):
 y=df.activity.to_numpy(float);p=df[pred].to_numpy(float);n=len(y);k=min(k,n);percent=(rankdata(y,method='average')-1)/(n-1) if n>1 else np.full(n,.5)
 cutoff=np.sort(p)[-k];selected=p>cutoff;tied=p==cutoff;remaining=k-int(selected.sum());mass=selected.astype(float)+tied*(remaining/tied.sum())
 return dict(n=n,k=k,top5_mean_percentile=float(mass@percent/k),selected_mean_activity=float(mass@y/k),random_tie_expectation=.5,cutoff_ties=int(tied.sum()),sequence_components=df.sequence_group.nunique())
