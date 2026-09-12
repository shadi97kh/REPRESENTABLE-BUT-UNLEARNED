"""Small deterministic correctness diagnostics; no sampling or fitted comparisons."""
import json
import math
from pathlib import Path
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from .core import (sinusoidal_profile, response_derivatives, reconstruct_profile,
                   preflight, transport_score, separation)

ACTIONS = [np.zeros(4), *np.eye(4)]
COEFFICIENTS = np.array([[.7,.8],[.9,.6]])
DIRECTION = np.array([[.1,-.08],[-.03,.04]])


def score_check():
    delta, step = .1, 1e-5
    profiles = [sinusoidal_profile(.04,1), sinusoidal_profile(-.03,1.3)]
    errors = []
    mean_scores = []
    def quantities(a,z):
        F,Fp,Fpp,x,y,hp = response_derivatives(a,z,delta,COEFFICIENTS,profiles)
        da,db = a[:2] @ DIRECTION
        r1,r2 = math.exp(x)*math.sin(x), .5*math.exp(y)*math.sin(1.3*y)
        r1p = math.exp(x)*(math.sin(x)+math.cos(x))
        r2p = .5*math.exp(y)*(math.sin(1.3*y)+1.3*math.cos(1.3*y))
        dF = r1+r2+profiles[0](x,1)*da+profiles[1](y,1)*db
        dFp = r1p+r2p*hp+profiles[0](x,2)*da+profiles[1](y,2)*db*hp
        return F,transport_score(z,Fp,Fpp,dF,dFp)
    for a in ACTIONS:
        for z in [-1.,.2,1.3]:
            y,score = quantities(a,z)
            logs=[]
            for t in [-step,step]:
                pt = [sinusoidal_profile(.04+t,1),sinusoidal_profile(-.03+.5*t,1.3)]
                ct = COEFFICIENTS+t*DIRECTION
                inverse = brentq(lambda zz:response_derivatives(a,zz,delta,ct,pt)[0]-y,-20,20,xtol=1e-13)
                Fp = response_derivatives(a,inverse,delta,ct,pt)[1]
                logs.append(-.5*inverse**2-.5*math.log(2*math.pi)-math.log(Fp))
            errors.append(abs(score-(logs[1]-logs[0])/(2*step)))
        val,err=quad(lambda z:quantities(a,z)[1]*math.exp(-z*z/2)/math.sqrt(2*math.pi),-12,12,epsabs=2e-10)
        mean_scores.append(dict(value=val,quadrature_error_estimate=err))
    return dict(fixed_outcome_derivative_max_abs_error=max(errors), number_of_points=len(errors),
                score_means=mean_scores,interval=[-12,12],status='NUMERICALLY CHECKED; not a DQM proof')


def reconstruction_check():
    profiles = [sinusoidal_profile(.04,1),sinusoidal_profile(-.03,1.3)]
    errors=[]
    for delta in [.1,.01]:
        curve=lambda a,z:response_derivatives(ACTIONS[a],z,delta,COEFFICIENTS,profiles)[0]
        for component in [1,2]:
            for x in [-1.2,.3,1.7]:
                for K in [2,5,8]:
                    recovered=reconstruct_profile(curve,component,x,delta,K)
                    expected=profiles[component-1](x)-profiles[component-1](x-K)
                    errors.append(abs(recovered-expected))
    return dict(telescoping_identity_max_abs_error=max(errors), number_of_values=len(errors),
                status='NUMERICALLY CHECKED using exact analytic quantile fixtures, no observations')


def near_null_check():
    # New fixed-delta unknown-profile directions, distinct from prior delta diagnostics.
    delta=.1
    profiles=[sinusoidal_profile(0,1),sinusoidal_profile(0,1)]
    output=[]
    for M in [2,4,6]:
        def relative(x,k):
            t=x+M
            polynomials=[1,1-2*t,4*t*t-4*t-1,-5+6*t+12*t*t-8*t**3]
            return M*(polynomials[k]*math.exp(-t*t)-math.exp(-M*M))
        norm2=0.
        norm_error=0.
        for left,right in [(-np.inf,-M),(-M,np.inf)]:
            val,err=quad(lambda x:sum(relative(x,k)**2 for k in range(4))/(1+x*x),left,right,epsabs=1e-9)
            norm2+=val;norm_error+=err
        norm=math.sqrt(norm2)
        info=0.; info_error=0.
        for a in ACTIONS:
            def integrand(z):
                _,Fp,Fpp,x,_,_=response_derivatives(a,z,delta,COEFFICIENTS,profiles)
                score=transport_score(z,Fp,Fpp,math.exp(x)*relative(x,0),math.exp(x)*relative(x,1))
                return score**2*math.exp(-z*z/2)/math.sqrt(2*math.pi)/5
            val,err=quad(integrand,-12,12,epsabs=1e-11,points=[-M,-M-1])
            info+=val;info_error+=err
        c1,c2=COEFFICIENTS[:,0]
        def derivative_integrand(z):
            r=lambda x:math.exp(x)*relative(x,0)
            return (r(z+c1+c2)-r(z+c1)-r(z+c2)+r(z))*math.exp(-z*z/2)/math.sqrt(2*math.pi)
        derivative,error=quad(derivative_integrand,-12,12,epsabs=1e-12,points=[-M,-M-1])
        output.append(dict(M=M,profile_norm=norm,score_norm_normalized=math.sqrt(info)/norm,
                           target_derivative_normalized=derivative/norm,
                           target_to_score_ratio=abs(derivative)/math.sqrt(info),
                           profile_norm2_quad_error=norm_error,score_norm2_quad_error=info_error,
                           target_derivative_quad_error=error))
    return dict(rows=output,allocation=[.2]*5,score_and_target_interval=[-12,12],
                interpretation='Both score and target derivatives shrink. This does not decide target regularity.',
                status='NUMERICALLY CHECKED; quadrature estimates are not certified bounds')


def run(out):
    out=Path(out);out.mkdir(parents=True,exist_ok=False)
    result=dict(reconstruction=reconstruction_check(),score=score_check(),near_null=near_null_check(),
                preflight=[dict(delta=d,**preflight(d,[1000]*5,2,2)) for d in [.1,.01]],
                sample_fits=0,random_draws=0,training_runs=0,gpu_use=0,
                note='Exact-function diagnostic fixtures are oracle evidence only, never learner inputs.')
    assert result['reconstruction']['telescoping_identity_max_abs_error']<1e-10
    assert result['score']['fixed_outcome_derivative_max_abs_error']<1e-6
    with (out/'diagnostics.json').open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
    print(json.dumps(result,indent=2,allow_nan=False))
