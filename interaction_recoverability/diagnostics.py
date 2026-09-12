"""Fixed deterministic quadrature checks only; no samples, fitting or benchmarks."""
import math
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from .model import s,h,coefficients,log_response,log_derivative,DELTAS,validate_delta
from .theory import experiment_bound,LIMIT

RADIUS=12.
EPSABS=2e-11
EPSREL=2e-11
PHI_CONSTANT=1/math.sqrt(2*math.pi)

def integrate(f,points=(0.,)):
    knots=sorted(set([-RADIUS,RADIUS]+[float(x) for x in points if -RADIUS<x<RADIUS]))
    parts=[quad(f,l,r,epsabs=EPSABS,epsrel=EPSREL,limit=150) for l,r in zip(knots[:-1],knots[1:])]
    return sum(v for v,e in parts),sum(e for v,e in parts)

def moment_tail(delta,power=1):
    """Analytic upper bound for omitted |Z|>12 integral of exp(power*h_delta)."""
    validate_delta(delta,boundary=True)
    if not np.isfinite(power) or power<=0:raise ValueError('positive moment order required')
    a=1-2*power*delta
    if a<=0:raise ValueError('requested exponential moment is not finite under this bound')
    positive=math.exp(power*delta/2+power**2/(2*a))/math.sqrt(a)*ndtr(-math.sqrt(a)*(RADIUS-power/a))
    negative=math.exp(power**2/2)*ndtr(-RADIUS-power)
    return float(positive+negative)

def moment(delta):
    value,error=integrate(lambda z:math.exp(float(h(delta,z,boundary=True))-.5*z*z)*PHI_CONSTANT)
    return dict(value=value,quadrature_error_estimate=error,analytic_tail_upper=moment_tail(delta))

def wasserstein(delta,x):
    validate_delta(delta)
    a,b=coefficients('T',x)-coefficients('Q',x)
    def f(z):
        # Preserve cancellation for the axis expression as delta tends to zero.
        difference=math.exp(z)*((a+b)+b*math.expm1(delta*float(s(z))))
        return abs(difference)*math.exp(-.5*z*z)*PHI_CONSTANT
    roots=[0.]
    if a*b<0:
        t=math.log(-a/b)/delta
        root=math.copysign(math.sqrt(2*t*t/(1+math.sqrt(1+4*t*t))),t)
        roots.append(root)
    value,error=integrate(f,roots)
    tail=abs(a)*moment_tail(0)+abs(b)*moment_tail(delta)
    return dict(value=value,quadrature_error_estimate=error,analytic_tail_upper=float(tail))

def hellinger(delta):
    """H^2 on axis e2 by density Jacobians; tails bounded in both laws."""
    x=np.array([0.,1.,0.,0.])
    def inverse_q(logy):
        return brentq(lambda w:float(log_response('Q',x,w,delta))-logy,-80.,80.,xtol=1e-12)
    def integrand(z):
        w=inverse_q(float(log_response('T',x,z,delta)))
        ratio_log=-.5*(w*w-z*z)+float(log_derivative('T',x,z,delta)-log_derivative('Q',x,w,delta))
        return math.expm1(.5*ratio_log)**2*math.exp(-.5*z*z)*PHI_CONSTANT
    value,error=integrate(integrand)
    wlo=inverse_q(float(log_response('T',x,-RADIUS,delta)))
    whi=inverse_q(float(log_response('T',x,RADIUS,delta)))
    tail=2*ndtr(-RADIUS)+ndtr(wlo)+ndtr(-whi)
    return dict(value=value,quadrature_error_estimate=error,analytic_tail_upper=float(tail),
                convention='H^2=int(sqrt(p)-sqrt(q))^2',inverse_xtol=1e-12,
                certified=False,limitation='Root/quadrature roundoff not enclosed; analytic tail formula evaluated in floating point')

def run():
    rows=[]
    for delta in DELTAS:
        m=moment(delta);axis=wasserstein(delta,[0,1,0,0]);joint=wasserstein(delta,[1,1,0,0])
        contrast_gap=(math.e-math.sqrt(math.e))*((math.e-1)*m['value']-(math.sqrt(math.e)-1)*math.sqrt(math.e))
        rows.append(dict(delta=delta,axis2_w1=axis,joint_w1=joint,ratio=joint['value']/axis['value'],
                         moment=m,primary_interaction_gap=contrast_gap,axis2_hellinger_squared=hellinger(delta),
                         complete_experiment_bound=experiment_bound(delta,[32]*5)))
    return dict(status='NUMERICALLY CHECKED',primary_target='I12=mu(e1+e2)-mu(e1)-mu(e2)+mu(0)',
                executed_deltas=list(DELTAS),integration_radius=RADIUS,epsabs=EPSABS,epsrel=EPSREL,
                limiting_gap=LIMIT,rows=rows,model_fits=0,training_steps=0,gpu_use=False,
                numerical_limit='Quadrature estimates and analytic tail bounds are reported separately; no certified numerical enclosure')
