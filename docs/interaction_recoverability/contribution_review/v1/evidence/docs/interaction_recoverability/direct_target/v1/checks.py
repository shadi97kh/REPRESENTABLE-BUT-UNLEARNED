"""Small deterministic analytic-identity checks; no samples, fitting or pilot imports."""
import json
import math
import os
from pathlib import Path
import sys
import time

ROOT=Path(__file__).resolve().parents[4]
sys.path.insert(0,str(ROOT))
import numpy as np
from scipy.integrate import quad
from interaction_recoverability import direct_target_v1 as d


def integrate(f, lower=-40., upper=90.):
    return quad(f,lower,upper,epsabs=2e-9,epsrel=2e-9,limit=180)[0]


def profile(x, component, oscillatory):
    frequency=1.3 if component==1 else 1.7
    return np.exp(x)*(1+(.025*np.sin(frequency*x) if oscillatory else 0.))


def main():
    start=time.process_time()
    assert len(os.sched_getaffinity(0))==1
    for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','NUMEXPR_NUM_THREADS']:
        assert os.environ.get(key)=='1',key
    assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
    assert not d.SAMPLED_DATA_FITTING_ENABLED
    results=[]
    for delta,oscillatory,c,K in [(.1,False,(.5,.75,.65,.9),1),(.03,True,(.6,.8,.7,.95),4)]:
        target=0.;direct=0.;remainder=0.;bound=0.
        for j,shifts in [(1,c[:2]),(2,c[2:])]:
            c1,c2=shifts
            g=lambda x:profile(x,j,oscillatory)
            latent=lambda z:z if j==1 else d.h(z,delta)
            tj=integrate(lambda z:float((g(latent(z)+c1+c2)-g(latent(z)+c1)-g(latent(z)+c2)+g(latent(z)))*d.phi(z)),-14,14)
            dj=integrate(lambda x:float(d.weight(x,c1,c2,j,delta,K)*(g(x+1)-g(x))))
            rj=integrate(lambda z:float((g(latent(z)+c1+c2-K)-g(latent(z)+c1-K)-g(latent(z)+c2-K)+g(latent(z)-K))*d.phi(z)),-14,14)
            moment=integrate(lambda z:float(np.exp(latent(z))*d.phi(z)),-14,14)
            bj=10*math.exp(-K)*moment*math.expm1(c1)*math.expm1(c2)
            assert abs(tj-dj-rj)<2e-7,(j,tj,dj,rj)
            assert abs(rj)<=bj+1e-10
            if not oscillatory:
                assert abs(tj-moment*math.expm1(c1)*math.expm1(c2))<1e-8
            target+=tj;direct+=dj;remainder+=rj;bound+=bj
        # Independent source-law expression, including reused q0 and h Jacobian.
        def source_expression(z):
            q0=profile(z,1,oscillatory)+profile(d.h(z,delta),2,oscillatory)
            q3=profile(z+1,1,oscillatory)+profile(d.h(z,delta),2,oscillatory)
            q4=profile(z,1,oscillatory)+profile(d.h(z,delta)+1,2,oscillatory)
            a0,a3,a4=d.latent_weights(z,c,delta,K)
            return float(a0*q0+a3*q3+a4*q4)
        source=integrate(source_expression,-30,16)
        assert abs(source-direct)<3e-7,(source,direct)
        grid=np.linspace(-6,6,21)
        w0,w3,w4=d.probability_weights_at_z(grid,c,delta,K,3)
        assert np.max(np.abs(w0+w3+w4))<1e-10
        assert all(np.all(w[np.abs(grid)>=4]==0) for w in [w0,w3,w4])
        inv=d.hinverse(d.h(grid,delta),delta)
        assert np.max(np.abs(inv-grid))<2e-14
        density_mass=integrate(lambda x:float(d.density(x,2,delta)),-90,90)
        assert abs(density_mass-1)<1e-8
        results.append(dict(delta=delta,oscillatory=oscillatory,coefficients=c,K=K,target=target,
            direct=direct,remainder=remainder,remainder_bound=bound,
            identity_absolute_error=abs(target-direct-remainder),source_weight_absolute_error=abs(source-direct),
            density_mass=density_mass,max_inverse_error=float(np.max(np.abs(inv-grid)))))
    # Periodic residue is tiny but nonzero. This is not an all-five-law lower bound.
    phase=.125
    measured=float(d.weight(-20+phase,.5,.5,1,.1,80))
    fourier=sum(8*math.exp(-2*math.pi**2*k*k)*math.cos(2*math.pi*k*phase) for k in range(1,8,2))
    assert abs(measured-fourier)<2e-15
    x=np.linspace(-4,4,17)
    actual=d.weight(x,1.,.7,1,.1,80)
    expected=d.phi(x-.7)-d.phi(x)
    integer_error=float(np.max(np.abs(actual-expected)))
    assert integer_error<2e-15
    # Derivative of the finite primitive, checked independently of density sum.
    eps=1e-5
    primitive_error=0.
    for component in [1,2]:
        fd=(d.weight(x+eps,.6,.8,component,.1,3,True)-d.weight(x-eps,.6,.8,component,.1,3,True))/(2*eps)
        primitive_error=max(primitive_error,float(np.max(np.abs(fd-d.weight(x,.6,.8,component,.1,3)))))
    assert primitive_error<2e-9
    record=dict(status='all deterministic checks passed',fixtures=results,
        periodic_residue=measured,fourier_residue=fourier,periodic_absolute_error=abs(measured-fourier),
        integer_shift_absolute_error=integer_error,primitive_derivative_error=primitive_error,
        measured_child_cpu_s=time.process_time()-start,affinity=sorted(os.sched_getaffinity(0)),
        workers=1,numerical_threads=1,gpu_s=0,
        limitations='Ordinary float64 diagnostics on stated finite integration domains, not certified enclosures. No sampled observations, fits, minimax tests or new pilot results.')
    with (Path(__file__).parent/'checks.json').open('x') as f:json.dump(record,f,indent=2)
    print(json.dumps(record,indent=2))


if __name__=='__main__':main()
