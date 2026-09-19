"""Frozen deterministic efficiency geometry; no sample input or fitting interface."""
import hashlib, importlib.util, json, math, os, resource, time
from pathlib import Path
import numpy as np
from scipy.special import ndtr
HERE=Path(__file__).resolve().parent
PRIOR=HERE.parents[1]/'practical_feasibility/v1'
spec=importlib.util.spec_from_file_location('saved_influence',PRIOR/'influence_variance.py')
m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
LD=np.longdouble

# Frozen before any output. Pair coefficient order, then k=1 and k=2 blocks.
DICTIONARY=['alpha1','beta1','alpha2','beta2']+[
    f'g{j}_{kind}_k{k}' for k in (1,2) for j in (1,2) for kind in ('sin','cos_minus_one')]
CONFIGS=[
    dict(name='coarse',depth=512,central_order=128,outer_step=.125,score_step=.25,score_order=8),
    dict(name='fine',depth=512,central_order=256,outer_step=.0625,score_step=.125,score_order=12)]
RADIUS=12
LATTICE=24
RANK_RELATIVE_TOLERANCE=1e-12
SOURCE_SHIFTS=np.asarray(m.SHIFTS,dtype=float)

def h(z):return z+.1*z*np.sqrt(1+z*z)
def hp(z):return 1+.1*(1+2*z*z)/np.sqrt(1+z*z)
def hpp(z):return .1*z*(3+2*z*z)/(1+z*z)**1.5

def profile(x,k,kind):
    ex=np.exp(x)
    if kind=='sin':return ex*np.sin(k*x),ex*(np.sin(k*x)+k*np.cos(k*x))
    return ex*(np.cos(k*x)-1),ex*(np.cos(k*x)-1-k*np.sin(k*x))

def score_geometry(cfg):
    zz,ww=m.gauss_grid(-LD(RADIUS),LD(RADIUS),LD(cfg['score_step']),cfg['score_order'])
    z=np.asarray(zz,dtype=float);weights=np.asarray(ww*m.phi(zz),dtype=float)
    hz=h(z);dh=hp(z);ddh=hpp(z)
    scores=[];sources=[]
    for a,(s,t) in enumerate(SOURCE_SHIFTS):
        x=z+s;y=hz+t;ex=np.exp(x);ey=np.exp(y)
        fp=ex+dh*ey;fpp=ex+(ddh+dh*dh)*ey
        q=np.zeros((len(z),12));dq=np.zeros_like(q)
        if a in (1,2):
            col=2*(a-1)
            q[:,col]=ex;dq[:,col]=ex
            q[:,col+1]=ey;dq[:,col+1]=dh*ey
        col=4
        for k in (1,2):
            for j in (1,2):
                for kind in ('sin','cos_minus_one'):
                    rr,dr=profile(x if j==1 else y,k,kind)
                    q[:,col]=rr;dq[:,col]=dr if j==1 else dh*dr
                    col+=1
        sc=z[:,None]*q/fp[:,None]-dq/fp[:,None]+q*(fpp/fp**2)[:,None]
        scores.append(sc);sources.append((sc.T*weights)@sc/5)
    G=sum(sources)
    # Direct four-shift target derivative, independent of the complex-moment check.
    b=np.zeros(12)
    al=np.asarray(m.ALPHA,dtype=float);be=np.asarray(m.BETA,dtype=float)
    m0=math.exp(.5);mh=float(weights@np.exp(hz))
    b[:4]=[m0*np.exp(al[0])*np.expm1(al[1]),mh*np.exp(be[0])*np.expm1(be[1]),
           m0*np.exp(al[1])*np.expm1(al[0]),mh*np.exp(be[1])*np.expm1(be[0])]
    col=4;complex_b=[]
    for k in (1,2):
        w=1+1j*k
        for j in (1,2):
            base=z if j==1 else hz;c=al if j==1 else be
            moment=np.exp(w*w/2) if j==1 else weights@np.exp(w*hz)
            dc=moment*np.expm1(w*c[0])*np.expm1(w*c[1])
            real_target=(m0 if j==1 else mh)*np.prod(np.expm1(c))
            for kind in ('sin','cos_minus_one'):
                diff=(profile(base+sum(c),k,kind)[0]-profile(base+c[0],k,kind)[0]
                      -profile(base+c[1],k,kind)[0]+profile(base,k,kind)[0])
                b[col]=weights@diff
                complex_b.append(float(dc.imag if kind=='sin' else dc.real-real_target))
                col+=1
    lower=[]
    for d in (4,8,12):
        gg=G[:d,:d];bb=b[:d]
        ev,U=np.linalg.eigh(gg)
        tol=RANK_RELATIVE_TOLERANCE*max(ev)
        retained=ev>tol
        # No ridge. Report all positive eigenmodes separately if rank tolerance drops any.
        assert min(ev)>0,'Numerical failure despite analytic positive definiteness'
        invcoef=np.linalg.solve(gg,bb)
        lower.append(dict(directions=d,bound=float(bb@invcoef),eigenvalues=ev.tolist(),
            minimum_eigenvalue=float(min(ev)),condition_number=float(max(ev)/min(ev)),
            rank_tolerance=float(tol),numerical_rank=int(sum(retained)),
            dropped_b_norm=float(np.linalg.norm((U.T@bb)[~retained])),
            thresholded_pseudoinverse_bound=float(np.sum((U.T@bb)[retained]**2/ev[retained])),
            solve_residual=float(np.linalg.norm(gg@invcoef-bb)),
            dual_coefficients=invcoef.tolist()))
    return dict(gram=G.tolist(),source_grams=[g.tolist() for g in sources],b=b.tolist(),
        lower_bounds=lower,score_mean_max=float(max(np.max(np.abs(weights@s)) for s in scores)),
        complex_target_derivative_max_error=float(np.max(np.abs(b[4:]-complex_b))),
        gram_tail_entry_bound=1e-28,
        tail_bound_scope='Analytic full-line tail bound at R=12; does not certify retained quadrature or roundoff.')

class MedianAccumulator(m.Atoms):
    """Stream linear source-zero covariance; never allocate sample or full atom bank."""
    def __init__(self):self.c=np.zeros(5,dtype=LD);self.nodes=0
    def add(self,a,z,w):
        if a!=0:return
        z=np.asarray(z,dtype=LD).reshape(-1);w=np.asarray(w,dtype=LD)
        if w.ndim==1:w=np.broadcast_to(w,(len(z),5))
        assert w.shape==(len(z),5)
        amp=-LD('.5')*m.fp(0,z)/m.phi(z)*np.asarray(ndtr(-np.abs(z).astype(float)),dtype=LD)
        self.c+=np.sum(w*amp[:,None],axis=0,dtype=LD);self.nodes+=len(z)

def median_geometry(cfg,prior):
    atom=MedianAccumulator();K=cfg['depth'];Q=cfg['central_order']
    # Split the known median kink and outward-weight join explicitly.
    for left,right in ((-LD(RADIUS),m.B),(m.B,LD(0)),(LD(0),LD(RADIUS))):
        z,wq=m.gauss_grid(left,right,LD(cfg['outer_step']))
        js=np.arange(1,LATTICE+1,dtype=LD) if left>=m.B else -np.arange(LATTICE,dtype=LD)
        W=np.sum(m.L(z[:,None]+js[None,:]),axis=1,dtype=LD)
        if left<m.B:W=-W
        v=np.zeros((len(z),5),dtype=LD);v[:,0]=wq*(m.A(z)-W);atom.add(0,z,v)
    t,wq=np.polynomial.legendre.leggauss(Q)
    t=m.B+(np.asarray(t,dtype=LD)+1)/2;wq=np.asarray(wq,dtype=LD)/2
    P=np.sum(m.L(t[:,None]+np.arange(-LATTICE,LATTICE+1,dtype=LD)),axis=1,dtype=LD)
    v=np.zeros((Q,5),dtype=LD);v[:,0]=P*wq
    atom.difference(t,-1,v,K)
    transform=np.eye(5,dtype=LD)
    for i in range(2):
        a,b=m.ALPHA[i],m.BETA[i]
        J=np.array([[np.exp(-1+a),np.exp(-m.h(LD(1))+b)],
                    [np.exp(1+a),np.exp(m.h(LD(1))+b)]],dtype=LD)
        det=J[0,0]*J[1,1]-J[0,1]*J[1,0]
        transform[1+2*i:3+2*i,1+2*i:3+2*i]=np.array([[J[1,1],-J[0,1]],[-J[1,0],J[0,0]]])/det
        for j,z in enumerate((LD(-1),LD(1))):
            unit=np.zeros(5,dtype=LD);unit[1+2*i+j]=1
            x=z+a;y=m.inv(m.h(z)+b)
            atom.add(i+1,[z],unit);atom.add(0,[y],-unit)
            ty=atom.transport(y,unit);tx=atom.transport(x,-unit)
            atom.difference([ty],[tx],unit,K)
    transformed=transform@atom.c
    gradient=np.asarray(prior['target_gradient_pair_order'],dtype=LD)
    c=transformed[0]+gradient@transformed[1:]
    V=LD(prior['complete_variance'])
    reduction=20*c*c;corrected=V-reduction
    assert 0<=reduction<V
    return dict(raw_anchor_and_residual_covariances=atom.c.astype(float).tolist(),
        anchor_and_coefficient_covariances=transformed.astype(float).tolist(),
        anchor_covariance=float(transformed[0]),coefficient_correction_covariance=float(gradient@transformed[1:]),
        complete_covariance=float(c),lambda_reference=float(20*c),
        reduction=float(reduction),reduction_percent=float(100*reduction/V),
        constructed_variance=float(V),corrected_variance=float(corrected),
        corrected_sd_constant=float(np.sqrt(corrected)),streamed_source0_nodes=atom.nodes,
        truncation_covariance_error_bound=5.905/math.sqrt(20),
        truncation_corrected_norm_error_bound=5.905,
        numeric_error_certified=False)

def main():
    output=HERE/'raw.json'
    if output.exists():raise FileExistsError('Retain every previous output')
    cpu=time.process_time();wall=time.monotonic()
    record=dict(status='started',dictionary=DICTIONARY,nested_sizes=[4,8,12],configs=CONFIGS,
        radius=RADIUS,lattice=LATTICE,pi=[.2]*5,normalization='N=5n; Hilbert weights 1/5',
        scope='Deterministic exact-reference geometry only; no samples/fits/simulations/backend.',
        numerical_rank_relative_tolerance=RANK_RELATIVE_TOLERANCE,results=[])
    try:
        assert len(os.sched_getaffinity(0))==1
        assert os.environ.get('OMP_NUM_THREADS')==os.environ.get('OPENBLAS_NUM_THREADS')=='1'
        assert os.environ.get('CUDA_VISIBLE_DEVICES')==''
        prior=json.loads((PRIOR/'variance_final.json').read_text())
        assert prior['status']=='passed'
        assert hashlib.sha256((PRIOR/'influence_variance.py').read_bytes()).hexdigest()==prior['code_sha256']
        record['prior_variance_sha256']=hashlib.sha256((PRIOR/'variance_final.json').read_bytes()).hexdigest()
        record['reused_influence_code_sha256']=prior['code_sha256']
        for cfg in CONFIGS:
            gg=score_geometry(cfg);mm=median_geometry(cfg,prior)
            record['results'].append(dict(config=cfg,score_geometry=gg,median=mm))
            assert gg['score_mean_max']<1e-9
            assert gg['complex_target_derivative_max_error']<1e-9
            bounds=[x['bound'] for x in gg['lower_bounds']]
            assert bounds==sorted(bounds) and bounds[-1]<mm['corrected_variance']
            print(json.dumps(dict(case=cfg['name'],lower_bounds=bounds,c=mm['complete_covariance'],
                                  corrected_variance=mm['corrected_variance'])),flush=True)
        coarse,fine=record['results']
        record['refinement_diagnostics']=dict(
            gram_operator_difference=float(np.linalg.norm(np.array(fine['score_geometry']['gram'])-coarse['score_geometry']['gram'],2)),
            b_l2_difference=float(np.linalg.norm(np.array(fine['score_geometry']['b'])-coarse['score_geometry']['b'])),
            covariance_difference=abs(fine['median']['complete_covariance']-coarse['median']['complete_covariance']),
            lower_bound_differences=[abs(y['bound']-x['bound']) for x,y in zip(coarse['score_geometry']['lower_bounds'],fine['score_geometry']['lower_bounds'])],
            interpretation='Refinement diagnostics are not numerical error certificates.')
        record['prior_gradient_max_difference']=float(np.max(np.abs(np.array(fine['score_geometry']['b'][:4])-prior['target_gradient_pair_order'])))
        record['gap']=dict(difference=fine['median']['corrected_variance']-fine['score_geometry']['lower_bounds'][-1]['bound'],
            ratio=fine['median']['corrected_variance']/fine['score_geometry']['lower_bounds'][-1]['bound'])
        record['status']='passed'
    except BaseException as ex:
        record.update(status='failed',error=repr(ex));raise
    finally:
        record.update(child_cpu_s=time.process_time()-cpu,child_wall_s=time.monotonic()-wall,
            peak_rss_kib=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,
            affinity=sorted(os.sched_getaffinity(0)),gpu_s=0,numerical_workers=1,numerical_threads=1,
            code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
        with output.open('x') as f:json.dump(record,f,indent=2)
if __name__=='__main__':main()

