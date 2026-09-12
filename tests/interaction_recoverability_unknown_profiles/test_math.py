"""New mathematical correctness checks; original diagnostic suite is not rerun."""
import math
import numpy as np
import pytest
from interaction_recoverability.unknown_profiles.core import (
    inverse_h,separation,preflight,error_bound,estimate_from_samples,
    sinusoidal_profile,response_derivatives,tail_envelope)
from interaction_recoverability.unknown_profiles.diagnostics import reconstruction_check,score_check
from interaction_recoverability.model import h


def test_private_anchor_telescoping_for_unknown_profiles():
    assert reconstruction_check()['telescoping_identity_max_abs_error']<1e-10


def test_density_score_differentiates_fixed_outcome_and_has_zero_mean():
    result=score_check()
    assert result['fixed_outcome_derivative_max_abs_error']<1e-6
    assert max(abs(x['value']) for x in result['score_means'])<1e-8


def test_warp_inverse_in_both_tails():
    for delta in [.01,.1]:
        for x in [-30,-1,0,1,30]:
            assert abs(float(h(delta,inverse_h(delta,x)))-x)<1e-10


def test_derived_global_coefficient_separation():
    profiles=[sinusoidal_profile(.04,1),sinusoidal_profile(-.03,1.3)]
    for delta in [.01,.1]:
        sep=separation(delta)
        for c in [np.array([[.51,.99],[.9,.6]]),np.array([[.99,.51],[.9,.6]]),np.array([[.6,.7],[.9,.6]])]:
            c0=np.array([[.75,.75],[.9,.6]])
            residual=[]
            for z,weight in [(-sep['R'],math.exp(sep['R'])),(sep['R'],math.exp(-sep['hR']))]:
                a=np.array([1,0,0,0])
                residual.append(weight*(response_derivatives(a,z,delta,c,profiles)[0]-response_derivatives(a,z,delta,c0,profiles)[0]))
            assert max(abs(x) for x in residual)>=sep['gamma']*np.max(np.abs(c[0]-c0[0]))-1e-10


def test_preflight_rejects_unobserved_quantile_precision():
    result=preflight(.1,[1000]*5,2,2)
    assert not result['informative']
    assert result['sufficient_log10_n']>30
    # Hypothetical counts evaluate an analytic bound; no samples are generated.
    result=error_bound(.1,[10**100]*5,2,2,20,20)
    assert result['informative']
    assert set(result['error_terms'])=={'stochastic_profile','profile_truncation','coefficient_and_grid','integration','outcome_tail'}


def test_no_sample_access_before_authorization_gate():
    class Trap:
        def __iter__(self):
            raise AssertionError('must not inspect any sample data')
    with pytest.raises(PermissionError):
        estimate_from_samples(Trap(),.1,K=2,B=2,grid_steps=2,integration_bins=2)


def test_boundary_and_invalid_inputs_fail_explicitly():
    for delta in [0,-.1,.2,float('nan')]:
        with pytest.raises(ValueError):separation(delta)
    with pytest.raises(ValueError):preflight(.1,[0]*5,2,2)
    with pytest.raises(ValueError):preflight(.1,[10]*5,True,2)
    assert tail_envelope(.1,8)<tail_envelope(.1,4)<tail_envelope(.1,0)
