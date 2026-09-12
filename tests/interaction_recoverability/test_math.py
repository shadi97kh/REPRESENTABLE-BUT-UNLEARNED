import math
import numpy as np
import pytest
from interaction_recoverability.model import response,h,s_prime,COEFFICIENTS,DELTAS
from interaction_recoverability.theory import LIMIT,experiment_bound,estimable_linear_functional,score_information_bound
from interaction_recoverability.estimator import (quantile_design,recover_coefficients_from_quantiles,
    primary_interaction,estimate_from_samples,quantile_error_bound)
from interaction_recoverability.diagnostics import moment_tail


def test_supports_monotonicity_and_normalization():
    for coefficients in COEFFICIENTS.values():
        assert [set(np.flatnonzero(row)) for row in coefficients]==[{0,1,2},{0,1,3}]
    z=np.array([-12.,-3.,0.,3.,12.])
    for delta in DELTAS:
        assert h(delta,0)==0
        assert np.all(1+delta*s_prime(z)>0)
        for label in ('T','Q'):assert np.all(np.diff(response(label,[1,1,0,0],z,delta))>0)


def test_axes_and_boundary_equalities():
    z=np.array([-3.,0.,2.])
    for i in range(4):
        x=np.eye(4)[i]*.37
        assert np.allclose(response('T',x,z,0,boundary=True),response('Q',x,z,0,boundary=True),rtol=1e-14)
        if i!=1:
            for delta in DELTAS:assert np.array_equal(response('T',x,z,delta),response('Q',x,z,delta))


def test_invalid_labels_actions_deltas_and_nonfinite_latents():
    for label in ('wrong','t',''):
        with pytest.raises(ValueError):response(label,[0]*4,0,.1)
    for x in ([0]*3,[0,0,0,2],[0,0,-1,0],[0,0,np.nan,0]):
        with pytest.raises(ValueError):response('T',x,0,.1)
    for d in (0,-.1,.1001,float('nan')):
        with pytest.raises(ValueError):response('T',[0]*4,0,d)
    with pytest.raises(ValueError):h(.1,float('inf'))


def test_limiting_expression_is_target_interaction_gap():
    mean=lambda label,x:float(response(label,x,0,0,boundary=True)*math.exp(.5))
    def interaction(label):return mean(label,[1,1,0,0])-mean(label,[1,0,0,0])-mean(label,[0,1,0,0])+mean(label,[0]*4)
    assert interaction('T')-interaction('Q')==pytest.approx(LIMIT,abs=1e-13)
    assert LIMIT==pytest.approx(1.88607083318,abs=5e-12)


def test_density_bound_full_experiment_counts():
    first=experiment_bound(.001,[32]*5)
    nuisance=experiment_bound(.001,[9999,9999,32,9999,9999])
    assert first==nuisance
    assert 0<first['tv_upper']<1 and first['mse_lower']>0
    assert experiment_bound(.001,[32,32,0,32,32])['h2_experiment_upper']==0
    assert 0<score_information_bound()<100
    with pytest.raises(ValueError):experiment_bound(.001,[32]*4)


@pytest.mark.parametrize('delta',DELTAS)
def test_exact_quantile_algebra_without_sample_fitting(delta):
    c=np.array([[math.exp(.5),math.e],[math.e,math.exp(.5)]])
    design=quantile_design(delta)
    q=c@design.T
    recovered=recover_coefficients_from_quantiles(q,delta)
    assert np.max(abs(recovered-c))<2e-10
    r=0.6744897501960817
    assert np.linalg.det(design)==pytest.approx(2*math.sinh(delta*r*math.sqrt(1+r*r)),rel=2e-11)
    assert primary_interaction(recovered,math.exp(.5))>0


def test_fitting_guard_and_finite_sample_bound_interface():
    with pytest.raises(PermissionError):estimate_from_samples([1,2],[3,4],.01,2.)
    assert not quantile_error_bound(.01,2)['informative']
    assert quantile_error_bound(.01,1000)['informative']


def test_estimability_target_avoids_nullspace_or_not():
    matrix=[[1.,1.],[0.,0.]]
    assert estimable_linear_functional(matrix,[1.,1.])['estimable']
    assert not estimable_linear_functional(matrix,[1.,-1.])['estimable']


def test_tail_bounds_and_invisible_misspecification():
    for d in DELTAS:
        assert 0<=moment_tail(d)<1e-18
        assert moment_tail(d,2)>0
    with pytest.raises(ValueError):moment_tail(.1,5)
    for x in np.vstack([np.zeros(4),np.eye(4)]):assert 3*x[0]*x[1]==0
    assert 3*np.array([1,1,0,0])[0]*np.array([1,1,0,0])[1]==3


def test_matching_lower_bound_and_protected_contrast_identity():
    from interaction_recoverability.theory import normalized_class_lower_bound
    from interaction_recoverability.estimator import protected_contrast_from_means
    small=normalized_class_lower_bound(.01,32)
    large=normalized_class_lower_bound(.01,32000)
    assert 0<large['mse_lower']<small['mse_lower']
    for model in ('T','Q'):
        # Pointwise identity implies the same expectation identity, with no fitting.
        z=np.array([-2.,0.,2.])
        base=response(model,[0,0,0,0],z,.01)
        e1=response(model,[1,0,0,0],z,.01)
        joint=response(model,[1,0,1,1],z,.01)
        anchors=response(model,[0,0,1,1],z,.01)
        assert np.allclose(joint-e1-anchors+base,(np.e-1)*(e1-base),rtol=1e-13)
    assert protected_contrast_from_means(3.,2.)==pytest.approx(np.e-1)


def test_recorded_quadrature_matches_report_with_declared_precision():
    import json
    from pathlib import Path
    p=Path('runs/interaction_recoverability/prep-20260912-082154/numerics-v1/diagnostics.json')
    data=json.loads(p.read_text())
    references=[(.556071,3.295503),(.042850,1.991465),(.004190,1.896344),(.000418,1.887096)]
    for row,(axis,joint) in zip(data['rows'],references):
        assert abs(row['axis2_w1']['value']-axis)<.500001e-6
        assert abs(row['joint_w1']['value']-joint)<.500001e-6
        assert row['axis2_hellinger_squared']['value']<=score_information_bound()*row['delta']**2
        assert row['primary_interaction_gap']>=LIMIT
