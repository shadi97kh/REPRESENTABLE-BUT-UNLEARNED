"""Fixed-reference deterministic checks of the analytic contraction construction.

Predeclared: this one reference; Neumann truncations 32,128,512; quadrature
orders 64,128 on fixed intervals; two analytic profile directions (trigonometric
and an independent Gaussian holdout), four coefficient directions and one mixed
direction. No fitted finite bank, samples, variance optimization or rank cutoff.
"""
import json
import math
import os
from pathlib import Path
import sys
import time
ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
import numpy as np
from numpy.polynomial.legendre import leggauss
from interaction_recoverability.direct_target_v1 import h, hprime, hinverse, phi, kernel

DELTA=.1
ALPHA=np.array([2/3,3/4]);BETA=np.array([3/5,4/5])
AA=np.array([0.,*ALPHA,1.,0.]);BB=np.array([0.,*BETA,0.,1.])
LAMBDA=575/613
LEVELS=(32,128,512)

def transform(x,j):return hinverse(h(x,DELTA)+j,DELTA)
def contract(x):return transform(x,2)-1
def cp(x):return hprime(x,DELTA)/hprime(transform(x,2),DELTA)

def profile(x,j,kind,derivative=False):
    x=np.asarray(x)
    if kind=='zero':return np.zeros_like(x)
    if kind=='trig':
        if j==1:
            v=np.sin(x);vp=np.cos(x)
        else:
            v=np.cos(1.3*x)-1;vp=-1.3*np.sin(1.3*x)
    else:
        center=.4 if j==1 else -.6
        width=1. if j==1 else .8
        eg=np.exp(-width*(x-center)**2)
        v=eg-math.exp(-width*center*center)
        vp=-2*width*(x-center)*eg
    return np.exp(x)*(v+vp if derivative else v)

def perturbation(a,z,kind,d,derivative=False):
    z=np.asarray(z);hz=h(z,DELTA);hp=hprime(z,DELTA)
    x=z+AA[a];y=hz+BB[a]
    value=profile(x,1,kind,derivative)+profile(y,2,kind,derivative)*(hp if derivative else 1)
    if a in (1,2):
        value=value+np.exp(x)*d[a-1]+np.exp(y)*d[a+1]*(hp if derivative else 1)
    return value

def evalue(x,kind,d,derivative=False):
    t1=transform(x,1);t2=transform(x,2);c=t2-1
    if not derivative:
        return perturbation(0,t2,kind,d)+perturbation(0,t1,kind,d)+perturbation(0,c,kind,d)-perturbation(4,x,kind,d)-perturbation(4,t1,kind,d)-perturbation(3,c,kind,d)
    hp=hprime(x,DELTA);t1p=hp/hprime(t1,DELTA);t2p=hp/hprime(t2,DELTA)
    return perturbation(0,t2,kind,d,True)*t2p+perturbation(0,t1,kind,d,True)*t1p+perturbation(0,c,kind,d,True)*t2p-perturbation(4,x,kind,d,True)-perturbation(4,t1,kind,d,True)*t1p-perturbation(3,c,kind,d,True)*t2p

def derivative_reconstruct(x,kind,d,levels=LEVELS):
    current=np.array(x,copy=True);product=np.ones_like(current);total=np.zeros_like(current);out={}
    for k in range(1,max(levels)+1):
        total-=product*evalue(current,kind,d,True)
        t=transform(current,2);product*=hprime(current,DELTA)/hprime(t,DELTA);current=t-1
        if k in levels:out[k]=total.copy()
    return out

def central_values(x,kind,d,terms=512):
    x=np.asarray(x);cur=np.concatenate([x.reshape(-1),[4.]])
    total=np.zeros(x.size)
    for k in range(terms):
        ev=evalue(cur,kind,d);total-=ev[:-1]-ev[-1]
        cur=contract(cur)
    base=sum(perturbation(3,float(j),kind,d)-perturbation(0,float(j),kind,d) for j in range(4))
    return (base+total).reshape(x.shape)

def reconstruct_values(points,kind,d):
    points=np.asarray(points);shifts=np.ceil(4-points).astype(int)
    central=points+shifts
    values=central_values(central,kind,d)
    for idx,(x,n) in enumerate(zip(points,shifts)):
        if n>=0:
            values[idx]-=sum(perturbation(3,x+j,kind,d)-perturbation(0,x+j,kind,d) for j in range(n))
        else:
            values[idx]+=sum(perturbation(3,x+n+j,kind,d)-perturbation(0,x+n+j,kind,d) for j in range(-n))
    return values

def nodes(order,left,right):
    x,w=leggauss(order)
    return (right+left)/2+(right-left)*x/2,(right-left)*w/2

def target_weights(z,primitive=False):
    return kernel(z,*ALPHA,1,DELTA,primitive)-((1 if primitive else hprime(z,DELTA))*kernel(h(z,DELTA),*BETA,2,DELTA,primitive))

def outward(z):
    z=np.asarray(z);positive=z>=4
    result=np.zeros_like(z)
    for k in range(24):
        result+=np.where(positive,target_weights(z+k+1),-target_weights(z-k))
    return result

def central_B(z):
    return -sum(target_weights(z+n,True)-target_weights(4.+n,True) for n in range(-24,25))

def check_profile_functional(kind,order):
    d=np.zeros(4);total=0.;direct=0.
    for left,right in [(-12.,4.),(4.,5.),(5.,12.)]:
        z,w=nodes(order,left,right)
        a=hprime(z,DELTA)*kernel(h(z,DELTA),*BETA,2,DELTA)
        ww=outward(z)
        total+=float(w@((a-ww)*perturbation(0,z,kind,d)+ww*perturbation(3,z,kind,d)))
        contrast=np.zeros_like(z)
        for j,base,c in [(1,z,ALPHA),(2,h(z,DELTA),BETA)]:
            contrast+=profile(base+c.sum(),j,kind)-profile(base+c[0],j,kind)-profile(base+c[1],j,kind)+profile(base,j,kind)
        direct+=float(w@(phi(z)*contrast))
    z,w=nodes(order,4.,5.);der=derivative_reconstruct(z,kind,d,(512,))[512]
    total+=float(w@(central_B(z)*der))
    return dict(order=order,source_functional=total,target_derivative=direct,absolute_error=abs(total-direct))

def main():
    start=time.process_time();assert len(os.sched_getaffinity(0))==1
    assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:assert os.environ.get(key)=='1'
    grid=np.linspace(3.,6.,33);cgrid=contract(grid)
    assert cgrid.min()>3 and cgrid.max()<6
    maxcp=float(cp(grid).max());assert maxcp<LAMBDA
    derivative_rows=[]
    for kind in ['trig','gaussian_holdout']:
        d=np.zeros(4)
        err=float(np.max(np.abs(evalue(grid,kind,d)-(profile(cgrid,1,kind)-profile(grid,1,kind)))))
        scale=1+max(float(np.max(np.abs(perturbation(a,zz,kind,d)))) for a,zz in [(0,transform(grid,2)),(0,transform(grid,1)),(0,cgrid),(4,grid),(4,transform(grid,1)),(3,cgrid)])
        print(json.dumps({'stage':'commutator','kind':kind,'absolute_error':err,'term_scale':scale,'scaled_error':err/scale}),flush=True)
        assert err/scale<2e-12,(err,scale)
        rec=derivative_reconstruct(grid,kind,d)
        errors={str(k):float(np.max(np.abs(rec[k]-profile(grid,1,kind,True)))) for k in LEVELS}
        assert errors['512']<2e-8,errors
        derivative_rows.append(dict(kind=kind,commutator_error=err,neumann_max_errors=errors))
    coeff_rows=[]
    for name,kind,d in [('trig_profile','trig',np.zeros(4)),('gaussian_holdout','gaussian_holdout',np.zeros(4))]+[(f'coefficient_{i}','zero',np.eye(4)[i]) for i in range(4)]+[('mixed','gaussian_holdout',np.array([.2,-.1,.15,-.05]))]:
        if kind=='zero':
            # Exact zero nuisance input makes its reconstruction exactly zero.
            recovered_profiles=np.zeros(8)
        else:
            queries=[]
            for i in range(2):
                queries.extend([-1+ALPHA[i],1+ALPHA[i],float(hinverse(h(-1.,DELTA)+BETA[i],DELTA)),float(hinverse(h(1.,DELTA)+BETA[i],DELTA))])
            recovered_profiles=reconstruct_values(np.array(queries),kind,d)
        recovered=np.zeros(4);determinants=[]
        for i in range(2):
            x1=recovered_profiles[4*i:4*i+2]
            z2=np.array([float(hinverse(h(z,DELTA)+BETA[i],DELTA)) for z in [-1.,1.]])
            r2=perturbation(0,z2,kind,d)-recovered_profiles[4*i+2:4*i+4]
            residual=perturbation(i+1,np.array([-1.,1.]),kind,d)-x1-r2
            matrix=np.array([[math.exp(z+ALPHA[i]),math.exp(float(h(z,DELTA))+BETA[i])] for z in [-1.,1.]])
            sol=np.linalg.solve(matrix,residual);recovered[i]=sol[0];recovered[i+2]=sol[1]
            determinants.append(float(np.linalg.det(matrix)))
        error=float(np.max(np.abs(recovered-d)));assert error<2e-7,(name,error)
        coeff_rows.append(dict(direction=name,expected=d.tolist(),recovered=recovered.tolist(),max_error=error,determinants=determinants))
    functional=[check_profile_functional(kind,order) for kind in ['trig','gaussian_holdout'] for order in [64,128]]
    assert max(r['absolute_error'] for r in functional)<2e-7,functional
    refinement=max(abs(functional[i]['source_functional']-functional[i+1]['source_functional']) for i in [0,2])
    out=dict(status='passed',reference=dict(delta=DELTA,alpha=ALPHA.tolist(),beta=BETA.tolist(),pi=[.2]*5),
        contraction=dict(sampled_min_C=float(cgrid.min()),sampled_max_C=float(cgrid.max()),sampled_max_derivative=maxcp,proved_upper_derivative=LAMBDA),
        derivative_checks=derivative_rows,coefficient_checks=coeff_rows,profile_functional_checks=functional,
        max_quadrature_refinement_difference=refinement,child_cpu_s=time.process_time()-start,
        numerical_workers=1,numerical_threads=1,affinity=sorted(os.sched_getaffinity(0)),gpu_s=0,
        scope='Analytic construction checks only. No finite score projection, no empirical samples, no learned model, no efficient variance estimate.',
        errors='Neumann tail controlled analytically in proof; lattice/response tail bounds separate there. 64/128 quadrature and float64 inverse roundoff are diagnostics, not certified integration enclosures.')
    with (Path(__file__).parent/'checks.json').open('x') as f:json.dump(out,f,indent=2)
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
