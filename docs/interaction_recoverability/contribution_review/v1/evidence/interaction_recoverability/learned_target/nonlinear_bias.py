"""Finite admissible paths and exact conditional bias identity diagnostics."""
import math
import numpy as np
import torch
from scipy.special import ndtr
from .profiles import K,ACTIONS,tensor,warp,relative_profile_directions,quadrature,target,inverse
from .geometry import score_bank

class PerturbedModel:
    def __init__(self,base,v,t):self.base=base;self.v=np.asarray(v);self.t=t;self.coef=base.coef.detach()+t*tensor(v[:4].reshape(2,2))
    def profile(self,j,x,order=0):
        weights=tensor(self.v[4+j*K:4+(j+1)*K])
        return self.base.profile(j,x,order)+self.t*torch.exp(x)*(relative_profile_directions(x)[order]@weights)
    def response(self,a,z,delta,derivatives=False):
        from .profiles import warp_derivatives
        aa=tensor(ACTIONS[a]);x=z+self.coef[0,0]*aa[0]+self.coef[1,0]*aa[1]+aa[2]
        y=warp(z,delta)+self.coef[0,1]*aa[0]+self.coef[1,1]*aa[1]+aa[3]
        F=self.profile(0,x)+self.profile(1,y)
        if not derivatives:return F
        hp,hpp=warp_derivatives(z,delta)
        return F,self.profile(0,x,1)+self.profile(1,y,1)*hp,self.profile(0,x,2)+self.profile(1,y,2)*hp*hp+self.profile(1,y,1)*hpp,x,y,hp

def admissible_length(model,v):
    coef=model.coef.detach().numpy().reshape(-1);v=np.asarray(v)
    nonzero=np.abs(v[:4])>1e-14
    tcoef=float(np.min(np.minimum(coef[nonzero]-.5,1-coef[nonzero])/np.abs(v[:4][nonzero]))) if nonzero.any() else float('inf')
    # Both source families have derivative ratios in [.6,1.4] and |g'''| <= 6.2.
    total=max(np.abs(v[4:4+K]).sum(),np.abs(v[4+K:]).sum())
    tprofile=min(.1/5,3.8/65)/total if total else float('inf')
    return .9*min(tcoef,tprofile)

def ks_statistic(model,samples,delta):
    values=[]
    for a,y in enumerate(samples):
        sorted_y=np.sort(y);p=ndtr(inverse(model,a,sorted_y,delta).detach().numpy());n=len(y)
        values.append(float(max(np.max(np.arange(1,n+1)/n-p),np.max(p-np.arange(n)/n))))
    return values

@torch.no_grad()
def finite_path(model,geo,c,v,delta,selection_samples=None):
    if v is None:return {'status':'no accepted enrichment direction','rows':[]}
    max_t=admissible_length(model,v);step=min(.05,.5*max_t)
    z,w=quadrature(128);w=w.numpy();I=float(target(model,delta,128))
    rows=[]
    for t in [-step,step]:
        alternative=PerturbedModel(model,v,t);means=[]
        for a in range(5):
            y=alternative.response(a,z,delta)
            psi=(score_bank(model,a,inverse(model,a,y,delta),delta)-geo['center'][a])@c
            means.append(float(w@psi))
        bias=float(target(alternative,delta,128))-I-float(geo['pi']@means)
        residual=float((geo['b']-geo['G']@c)@v)
        row=dict(t=t,admissible_limit=max_t,B=bias,linear_prediction=t*residual,
                 nonlinear_remainder=bias-t*residual,status='finite-path witness; not a class supremum upper bound')
        if selection_samples is not None:
            ks=ks_statistic(alternative,selection_samples,delta)
            eps=[math.sqrt(math.log(10/.05)/(2*len(y))) for y in selection_samples]
            row['source_compatible']=all(k<=e for k,e in zip(ks,eps));row['ks']=ks
        rows.append(row)
    return {'rows':rows,'status':'admissible coupled profile/coefficient path diagnostics'}

def conditional_mse(bias,stratum_variances,counts,weights):
    return float(bias*bias+sum(p*p*v/n for p,v,n in zip(weights,stratum_variances,counts)))
