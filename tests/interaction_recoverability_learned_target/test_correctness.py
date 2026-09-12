import ast
import math
from pathlib import Path
import numpy as np
import torch
from interaction_recoverability.learned_target.profiles import ProfileModel,tensor,target,negative_log_likelihood,inverse,quadrature
from interaction_recoverability.learned_target.geometry import DIM,score_bank,target_gradient,geometry,metric
from interaction_recoverability.learned_target.correction import regularized_solve,evaluation_correction,fixed_space
from interaction_recoverability.learned_target.nonlinear_bias import PerturbedModel,admissible_length,conditional_mse
from interaction_recoverability.learned_target.experiment import stratified_roles,load_source

def test_fixed_outcome_likelihood_gradient_includes_inverse_and_jacobian():
    model=ProfileModel('basis');delta=.1;z=tensor([-.7,.3,1.2]);a=5
    with torch.no_grad():model.theta[0,1]=.3;model.theta[1,9]=-.2
    y=model.response(a,z,delta).detach();loss=negative_log_likelihood(model,a,y,delta).mean();loss.backward()
    score=score_bank(model,a,z,delta).mean(0)
    assert np.max(np.abs(model.coef.grad.numpy().reshape(-1)+score[:4]))<1e-9
    assert np.max(np.abs(model.theta.grad.numpy().reshape(-1)+.08*score[4:]))<1e-9
    original=float(model.coef[0,0]);vals=[]
    for step in [-1e-5,1e-5]:
        with torch.no_grad():model.coef[0,0]=original+step
        vals.append(float(negative_log_likelihood(model,a,y,delta).mean()))
    assert abs((vals[1]-vals[0])/2e-5-float(model.coef.grad[0,0]))<1e-7

def test_density_normalization_and_mean_zero_scores():
    model=ProfileModel('neural',seed=7);delta=.03
    t,w=np.polynomial.legendre.leggauss(256);t=21*t+6;w=21*w
    y=tensor(np.exp(t));z=inverse(model,2,y,delta)
    density=torch.exp(-negative_log_likelihood(model,2,y,delta)).detach().numpy()
    assert abs(float(w@(density*np.exp(t)))-1)<1e-7
    zz,ww=quadrature(192)
    for a in range(5):assert np.max(np.abs(ww.numpy()@score_bank(model,a,zz,delta)))<1e-7

def test_profile_constraints_and_normalization_after_projection():
    for kind in ['neural','basis']:
        model=ProfileModel(kind,seed=11)
        with torch.no_grad():
            for p in model.parameters():p.mul_(15)
            model.project()
        x=tensor(np.linspace(-12,12,151))
        for j in range(2):
            assert abs(float(model.profile(j,tensor(0.)))-1)<1e-14
            ratio=(model.profile(j,x,1)/torch.exp(x)).detach().numpy()
            assert ratio.min()>=.5 and ratio.max()<=2
            for d in [2,3]:assert float((model.profile(j,x,d)/torch.exp(x)).abs().max())<10
    # These grid checks guard code; the global proof is algebraic in method.md.

def test_all_target_derivatives_and_coupled_finite_path():
    model=ProfileModel('neural',seed=5);delta=.1;b=target_gradient(model,delta,192)
    for k in [0,1,2,3,5,4+16+9]:
        v=np.zeros(DIM);v[k]=1;step=min(1e-5,admissible_length(model,v)/10)
        diff=(float(target(PerturbedModel(model,v,step),delta,192))-float(target(PerturbedModel(model,v,-step),delta,192)))/(2*step)
        assert abs(diff-b[k])<1e-7
    v=np.zeros(DIM);v[:4]=[.1,-.07,.05,.1];v[5]=.2;v[4+16+9]=-.15
    M=np.eye(DIM);geo=geometry(model,delta,np.ones(5)/5,M,192)
    c=regularized_solve(geo,fixed_space(),.01)['c'];z,w=quadrature(192);w=w.numpy()
    I=float(target(model,delta,192));Bs=[]
    for step in [-1e-4,1e-4,.01]:
        alt=PerturbedModel(model,v,step);Epsi=0
        for a in range(5):
            yy=alt.response(a,z,delta).detach()
            values=(score_bank(model,a,inverse(model,a,yy,delta),delta)-geo['center'][a])@c
            Epsi+=float(w@values)/5
        Bs.append(float(target(alt,delta,192))-I-Epsi)
    residual=float((geo['b']-geo['G']@c)@v)
    assert abs((Bs[1]-Bs[0])/2e-4-residual)<1e-6
    assert abs(Bs[2]-.01*residual)>1e-9

def test_linear_inverse_risk_and_unidentified_target():
    G=np.diag([1.,.03**2]);b=np.array([1.,2.]);lam=.01;n=101;r=.2
    c=np.linalg.solve(G+lam*np.eye(2),b);res=b-G@c
    theta=r*res/np.linalg.norm(res)
    exact=(res@theta)**2+c@G@c/n
    spectral=np.sum(np.diag(G)*b*b/(np.diag(G)+lam)**2)/n+r*r*np.sum(lam*lam*b*b/(np.diag(G)+lam)**2)
    assert abs(exact-spectral)<1e-14
    nullgeo={'G':np.diag([1.,0.]),'b':np.ones(2),'M':np.eye(2)}
    result=regularized_solve(nullgeo,np.eye(2),.1)
    assert result['score_null_target_norm']==1 and result['residual'][1]==1

def test_coupled_cancellation_does_not_hide_unidentified_target():
    model=ProfileModel('basis')
    with torch.no_grad():model.coef.copy_(tensor([[.75,.75],[.65,.85]]))
    z,w=quadrature(96);v=np.zeros(DIM);v[:2]=[1,-1]
    for a in range(5):assert np.max(np.abs(score_bank(model,a,z,0.)@v))<1e-13
    assert abs(target_gradient(model,0.)@v)>.1

def test_unequal_allocation_variance_and_conditional_bias_sign():
    counts=[7,11,5,13,17];scores=[np.arange(n,dtype=float)[:,None] for n in counts];weights=np.array([.1,.2,.3,.15,.25])
    correction,var=evaluation_correction(scores,np.ones(1),weights)
    assert abs(correction-sum(p*(n-1)/2 for p,n in zip(weights,counts)))<1e-12
    variances=[np.var(s[:,0],ddof=1) for s in scores]
    assert conditional_mse(-.3,variances,counts,weights)==.09+var

def test_roles_truth_interfaces_and_witness_algebra(tmp_path):
    samples=[np.arange(n) for n in [31,32,33,34,35]];parts=stratified_roles(samples,7)
    for n,p in zip(map(len,samples),parts):
        assert set(p[0]).isdisjoint(p[1]) and set(p[1]).isdisjoint(p[2]) and set(p[0]).isdisjoint(p[2])
        assert sorted(sum(p,[]))==list(range(n))
    np.savez(tmp_path/'bad.npz',y0=np.ones(10),truth=np.ones(4))
    try:load_source(tmp_path/'bad.npz');assert False
    except ValueError:pass
    for name in ['profiles.py','source_fit.py','geometry.py','correction.py','enrichment.py','nonlinear_bias.py','baselines.py','experiment.py']:
        tree=ast.parse((Path('interaction_recoverability/learned_target')/name).read_text())
        assert not any(isinstance(node,ast.ImportFrom) and node.module and 'evaluate' in node.module for node in ast.walk(tree))
    A=np.array([[2.,.2],[.2,1.]]);R=np.array([1.,3.]);v=np.linalg.solve(A,R);v/=np.linalg.norm(v)
    assert abs((R@v)**2/(v@A@v)-R@np.linalg.solve(A,R))<1e-12
