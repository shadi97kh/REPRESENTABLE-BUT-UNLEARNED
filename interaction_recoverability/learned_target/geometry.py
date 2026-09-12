"""Actual fixed-outcome score geometry in an independent analytic direction sieve."""
import math
import numpy as np
import torch
from scipy.special import roots_legendre
from .profiles import K,ACTIONS,relative_profile_directions,tensor,quadrature,warp

DIM=4+2*K
# Four coefficient directions and two initial profile shapes per component.
INITIAL=[0,1,2,3,4+9,4+14,4+K+9,4+K+14]

@torch.no_grad()
def score_bank(model,a,z,delta):
    z=tensor(z);_,Fp,Fpp,x,y,hp=model.response(a,z,delta,True)
    r1=relative_profile_directions(x);r2=relative_profile_directions(y)
    rows=[];primes=[];aa=ACTIONS[a]
    for i,j in [(0,0),(0,1),(1,0),(1,1)]:
        xx=x if j==0 else y;mult=1 if j==0 else hp
        rows.append(aa[i]*model.profile(j,xx,1));primes.append(aa[i]*model.profile(j,xx,2)*mult)
    dotF=torch.cat([torch.stack(rows,-1),torch.exp(x)[:,None]*r1[0],torch.exp(y)[:,None]*r2[0]],-1)
    dotFp=torch.cat([torch.stack(primes,-1),torch.exp(x)[:,None]*r1[1],(torch.exp(y)*hp)[:,None]*r2[1]],-1)
    return (z[:,None]*dotF/Fp[:,None]-dotFp/Fp[:,None]+dotF*(Fpp/Fp**2)[:,None]).numpy()

@torch.no_grad()
def target_gradient(model,delta,nodes=96):
    z,w=quadrature(nodes);w=w.numpy();b=np.zeros(DIM)
    for a,sign in [(5,1),(1,-1),(2,-1),(0,1)]:
        _,_,_,x,y,_=model.response(a,z,delta,True);aa=ACTIONS[a]
        for k,(i,j) in enumerate([(0,0),(0,1),(1,0),(1,1)]):
            b[k]+=sign*float(tensor(w)@(aa[i]*model.profile(j,x if j==0 else y,1)))
        b[4:4+K]+=sign*w@(torch.exp(x)[:,None]*relative_profile_directions(x)[0]).numpy()
        b[4+K:]+=sign*w@(torch.exp(y)[:,None]*relative_profile_directions(y)[0]).numpy()
    return b

def metric(nodes=4096):
    """Numerical original weighted Sobolev Gram, x=tan(t); errors remain explicit."""
    theta,w=roots_legendre(nodes);theta*=math.pi/2;w*=math.pi/2
    with torch.no_grad():r=[x.numpy() for x in relative_profile_directions(tensor(np.tan(theta)))]
    block=sum(v.T@(w[:,None]*v) for v in r)
    M=np.eye(DIM);M[4:4+K,4:4+K]=block;M[4+K:,4+K:]=block
    return (M+M.T)/2

def whitening(M,tolerance=1e-9):
    d,U=np.linalg.eigh((M+M.T)/2)
    keep=d>tolerance*max(1.,float(d.max()))
    if not keep.all():raise ArithmeticError('direction metric redundant: revise declared basis, do not silently erase target')
    return U@np.diag(1/np.sqrt(d))@U.T

def geometry(model,delta,pi,M,nodes=96):
    z,w=quadrature(nodes);w=w.numpy();raw=np.stack([score_bank(model,a,z,delta) for a in range(5)])
    center=np.einsum('j,ajk->ak',w,raw);S=raw-center[:,None,:]
    G=np.einsum('a,j,ajk,ajl->kl',pi,w,S,S)
    return {'G':(G+G.T)/2,'b':target_gradient(model,delta,nodes),'M':M,
            'center':center,'nodes':nodes,'pi':np.asarray(pi)}

def fit_geometry_error(model,delta,pi,M,nodes=96):
    lo=geometry(model,delta,pi,M,nodes);hi=geometry(model,delta,pi,M,2*nodes)
    return lo,{'gram_difference_norm':float(np.linalg.norm(lo['G']-hi['G'],2)),
               'gradient_difference_norm':float(np.linalg.norm(lo['b']-hi['b'])),
               'centering_max_abs':float(np.abs(lo['center']).max()),
               'status':'quadrature refinement differences, not certified upper errors'}
