"""Deterministic reference influence integration. No samples or estimator fits.

Sorted centered step kernels integrate the Brownian-bridge covariance without
an atom-by-atom covariance matrix. Extended precision is diagnostic, not an
interval certificate. Every configuration writes a new result exclusively.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import resource
import time
import numpy as np
from scipy.special import ndtr
from scipy.integrate import quad

LD=np.longdouble
ROOT=Path(__file__).resolve().parent
B=LD(-7)/4
ALPHA=np.array([LD(2)/3,LD(3)/4])
BETA=np.array([LD(3)/5,LD(4)/5])
SHIFTS=[(0,0),(ALPHA[0],BETA[0]),(ALPHA[1],BETA[1]),(1,0),(0,1)]
INV_RESIDUAL=LD(0)


def h(x):
    x=np.asarray(x,dtype=LD)
    return x+x*np.sqrt(1+x*x)/10


def hp(x):
    return 1+(1+2*x*x)/(10*np.sqrt(1+x*x))


def inv(y):
    global INV_RESIDUAL
    y=np.asarray(y,dtype=LD)
    x=2*y/(1+np.sqrt(1+LD('.4')*np.abs(y)))
    for _ in range(12):
        x=x-(h(x)-y)/hp(x)
    INV_RESIDUAL=max(INV_RESIDUAL,np.max(np.abs(h(x)-y),initial=LD(0)))
    return x


def phi(x):
    return np.exp(-x*x/2)/np.sqrt(2*np.arccos(LD(-1)))


def fp(a,z):
    s,t=SHIFTS[a]
    return np.exp(z+s)+hp(z)*np.exp(h(z)+t)


def signed_density(x,a,b,component):
    def density(t):
        if component==1:
            return phi(t)
        z=inv(t)
        return phi(z)/hp(z)
    return density(x-a-b)-density(x-a)-density(x-b)+density(x)


def A(x):
    return hp(x)*signed_density(h(x),*BETA,2)


def L(x):
    return signed_density(x,*ALPHA,1)-A(x)


def gauss_grid(left,right,step,order=4):
    count=math.ceil(float((right-left)/step))
    knots=np.linspace(LD(left),LD(right),count+1,dtype=LD)
    x,w=np.polynomial.legendre.leggauss(order)
    z=(knots[:-1,None]+knots[1:,None])/2+(knots[1:]-knots[:-1])[:,None]*np.array(x,dtype=LD)/2
    ww=np.broadcast_to((knots[1:]-knots[:-1])[:,None]*np.array(w,dtype=LD)/2,z.shape)
    return z.ravel(),ww.ravel()


class Atoms:
    def __init__(self):
        self.z=[[] for _ in range(5)]
        self.w=[[] for _ in range(5)]

    def add(self,a,z,w):
        z=np.asarray(z,dtype=LD).reshape(-1)
        w=np.asarray(w,dtype=LD)
        if w.ndim==1:
            w=np.broadcast_to(w,(len(z),5))
        assert w.shape==(len(z),5)
        self.z[a].append(z.copy()); self.w[a].append(w.copy())

    def E(self,x,w):
        ts=[x]+[inv(h(x)+j) for j in range(1,6)]
        for j in range(1,6):self.add(0,ts[j],w)
        for j in range(5):self.add(4,ts[j],-w)
        for k in range(1,5):
            self.add(0,ts[5]-k,w); self.add(3,ts[5]-k,-w)

    def difference(self,x,y,w,K):
        x=np.asarray(x,dtype=LD).reshape(-1)
        y=np.broadcast_to(np.asarray(y,dtype=LD),x.shape).copy()
        for _ in range(K):
            self.E(x,-w); self.E(y,w)
            x,y=inv(h(x)+5)-4,inv(h(y)+5)-4

    def transport(self,x,w):
        m=math.ceil(float(B-x))
        if m<0:
            for k in range(1,1-m):self.add(3,[x-k],w);self.add(0,[x-k],-w)
        if m>0:
            for k in range(m):self.add(3,[x+k],-w);self.add(0,[x+k],w)
        return x+m


def probability(z):
    z64=np.asarray(z,dtype=float)
    return np.asarray(ndtr(z64),dtype=LD),np.asarray(ndtr(-z64),dtype=LD)


def interval_masses(z):
    p,sf=probability(z)
    dz=np.diff(z); left=z[:-1];right=z[1:]
    mass=np.empty(len(z)+1,dtype=LD);mass[0]=p[0];mass[-1]=sf[-1]
    interior=np.where(left>=0,sf[:-1]-sf[1:],p[1:]-p[:-1])
    small=dz<LD('.001')
    # Five-node Gauss avoids subtracting almost equal CDF probabilities.
    x,w=np.polynomial.legendre.leggauss(5)
    mid=(left[small]+right[small])/2
    half=dz[small]/2
    interior[small]=half*np.sum(phi(mid[:,None]+half[:,None]*x)*w,axis=1,dtype=LD)
    mass[1:-1]=interior
    assert np.min(mass)>=0
    return p,sf,mass


def step_covariance(z,w,precision='extended'):
    idx=np.argsort(z,kind='stable');z=z[idx];w=w[idx]
    # Merge equal query locations before sums to retain exact algebraic pairing.
    start=np.r_[True,z[1:]!=z[:-1]]
    indices=np.flatnonzero(start)
    w=np.add.reduceat(w,indices,axis=0);z=z[indices]
    p,sf,mass=interval_masses(z)
    typ=LD if precision=='extended' else np.float64
    w=w.astype(typ);p=p.astype(typ);sf=sf.astype(typ);mass=mass.astype(typ)
    prefix=np.vstack([np.zeros((1,5),dtype=typ),np.cumsum(w*p[:,None],axis=0,dtype=typ)])
    suffix=np.vstack([np.cumsum((-w*sf[:,None])[::-1],axis=0,dtype=typ)[::-1],np.zeros((1,5),dtype=typ)])
    values=prefix+suffix
    cov=np.empty((5,5),dtype=typ)
    for j in range(5):
        for k in range(j,5):
            cov[j,k]=cov[k,j]=5*np.sum(mass*values[:,j]*values[:,k],dtype=typ)
    mean=np.sum(mass[:,None]*values,axis=0,dtype=typ)
    return cov,{'distinct_nodes':len(z),'mass_error':float(abs(np.sum(mass)-1)),
                'mean_kernel':np.asarray(mean,dtype=float).tolist()}


def main(case,K,Q,R,step,precision):
    path=ROOT/(case+'.json')
    if path.exists():raise FileExistsError('Previous diagnostic must be retained')
    wall=time.monotonic();cpu=time.process_time()
    record={'status':'started','case':case,'depth':K,'central_order':Q,'response_radius':R,
            'outer_panel_width':step,'lattice':24,'accumulator':precision,
            'numpy_longdouble_bits':np.finfo(LD).nmant+1,
            'scope':'Reference influence approximation, no sampled data. Refinements are not interval certificates.'}
    try:
        assert 16<=K<=512 and 16<=Q<=512 and 4<=R<=12 and .02<=step<=.5
        # Preflight size: at most 36 K Q + 144 K + O(R/step) atoms.
        cap=36*K*Q+144*K+10000
        assert cap<10000000
        record['preflight_atom_bound']=cap
        record['preflight_storage_upper_bytes']=cap*256
        atom=Atoms()
        for left,right in ((-LD(R),B),(B,LD(R))):
            z,wq=gauss_grid(left,right,LD(step))
            js=np.arange(1,25,dtype=LD) if left>=B else -np.arange(24,dtype=LD)
            W=np.sum(L(z[:,None]+js[None,:]),axis=1,dtype=LD)
            if left<B:W=-W
            v=np.zeros((len(z),5),dtype=LD);v[:,0]=wq*(A(z)-W);atom.add(0,z,v)
            v=np.zeros((len(z),5),dtype=LD);v[:,0]=wq*W;atom.add(3,z,v)
        t,wq=np.polynomial.legendre.leggauss(Q)
        t=B+(np.asarray(t,dtype=LD)+1)/2;wq=np.asarray(wq,dtype=LD)/2
        P=np.sum(L(t[:,None]+np.arange(-24,25,dtype=LD)),axis=1,dtype=LD)
        v=np.zeros((Q,5),dtype=LD);v[:,0]=P*wq
        atom.difference(t,-1,v,K)
        record['finite_periodization_integral_approx']=float(np.sum(P*wq,dtype=LD))
        Js=[]
        for i in range(2):
            a,b=ALPHA[i],BETA[i]
            Js.append(np.array([[np.exp(-1+a),np.exp(-h(LD(1))+b)],
                                [np.exp(1+a),np.exp(h(LD(1))+b)]],dtype=LD))
            for j,z in enumerate((LD(-1),LD(1))):
                unit=np.zeros(5,dtype=LD);unit[1+2*i+j]=1
                x=z+a;y=inv(h(z)+b)
                atom.add(i+1,[z],unit);atom.add(0,[y],-unit)
                # ζ=κ_i−κ_0(y)+u_y−u_x; normalization and common reference cancel.
                ty=atom.transport(y,unit);tx=atom.transport(x,-unit)
                atom.difference([ty],[tx],unit,K)
        transform=np.eye(5,dtype=LD)
        for i,J in enumerate(Js):
            det=J[0,0]*J[1,1]-J[0,1]*J[1,0]
            transform[1+2*i:3+2*i,1+2*i:3+2*i]=np.array([[J[1,1],-J[0,1]],[-J[1,0],J[0,0]]])/det
        mh,quaderr=quad(lambda z:math.exp(float(h(LD(z)))-z*z/2)/math.sqrt(2*math.pi),-12,12,epsabs=1e-11)
        m0=np.exp(LD('.5'))
        target=m0*np.prod(np.expm1(ALPHA))+mh*np.prod(np.expm1(BETA))
        bvec=np.array([m0*np.exp(ALPHA[0])*np.expm1(ALPHA[1]),
                       mh*np.exp(BETA[0])*np.expm1(BETA[1]),
                       m0*np.exp(ALPHA[1])*np.expm1(ALPHA[0]),
                       mh*np.exp(BETA[1])*np.expm1(BETA[0])],dtype=LD)
        complete=np.r_[LD(1),bvec]
        cov_total=np.zeros((5,5),dtype=LD);source=[];raw=[]
        profile_response=np.zeros(5,dtype=LD)
        for a in range(5):
            z=np.concatenate(atom.z[a]);w=np.concatenate(atom.w[a]);raw.append(len(z))
            sa,sb=SHIFTS[a]
            q=np.exp(z+sa)*np.sin(z+sa)+np.exp(h(z)+sb)*np.sin(2*(h(z)+sb))
            profile_response+=np.sum(w*q[:,None],axis=0,dtype=LD)
            amp=w*(fp(a,z)/phi(z))[:,None]
            cv,diagnostic=step_covariance(z,amp,precision)
            cv=transform@cv@transform.T;cov_total+=cv
            source.append({'source':a,'covariance_anchor_and_coefficients':cv.astype(float).tolist(),
                'target_variance_contribution':float(complete@cv@complete),
                'anchor_variance':float(cv[0,0]),'correction_variance':float(bvec@cv[1:,1:]@bvec),
                'anchor_correction_covariance':float(cv[0,1:]@bvec),**diagnostic})
        variance=complete@cov_total@complete
        record.update(status='passed',target_reference=float(target),target_gradient_pair_order=bvec.astype(float).tolist(),
            covariance_order=['anchor','alpha1','beta1','alpha2','beta2'],
            complete_variance=float(variance),asymptotic_sd_constant=float(np.sqrt(variance)),
            coefficient_covariance=cov_total[1:,1:].astype(float).tolist(),
            anchor_variance=float(cov_total[0,0]),correction_variance=float(bvec@cov_total[1:,1:]@bvec),
            anchor_correction_covariance=float(cov_total[0,1:]@bvec),sources=source,
            asymptotic_se={str(N):float(np.sqrt(variance/N)) for N in (1000,10000,100000,1000000)},
            profile_tangent_coefficient_residual=(transform@profile_response)[1:].astype(float).tolist(),
            inverse_max_equation_residual=float(INV_RESIDUAL),reference_mean_quad_error_estimate=quaderr,
            atoms_by_source=raw)
        # A small independent Brownian-bridge covariance identity check.
        zz=np.array([-1,LD('.2'),1],dtype=LD);ww=np.zeros((3,5),dtype=LD)
        ww[:,0]=[2,-3,4];cv,_=step_covariance(zz,ww)
        pp=ndtr(np.asarray(zz,dtype=float));direct=5*np.array([2,-3,4])@(np.minimum.outer(pp,pp)-np.outer(pp,pp))@np.array([2,-3,4])
        record['small_bridge_formula_error']=float(abs(cv[0,0]-direct))
        assert record['small_bridge_formula_error']<1e-11
        assert np.max(np.abs(record['profile_tangent_coefficient_residual']))<1e-8
        assert variance>0
    except BaseException as ex:
        record.update(status='failed',error=repr(ex));raise
    finally:
        record.update(child_cpu_s=time.process_time()-cpu,child_wall_s=time.monotonic()-wall,
                      child_peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
                      code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        path.write_text(json.dumps(record,indent=2)+'\n')
        print(json.dumps({k:record[k] for k in ('case','status','complete_variance','asymptotic_sd_constant','child_cpu_s') if k in record}))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('case');p.add_argument('--depth',type=int,default=128)
    p.add_argument('--order',type=int,default=64);p.add_argument('--radius',type=float,default=8)
    p.add_argument('--step',type=float,default=.25);p.add_argument('--precision',choices=['extended','float64'],default='extended')
    a=p.parse_args();main(a.case,a.depth,a.order,a.radius,a.step,a.precision)
