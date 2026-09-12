import copy
from dataclasses import replace
import numpy as np
import pytest
import torch
from intervention.transforms import (molecular_input,replace_mapped_atom,ProductLibrary,
    FrozenQueryOracle,rebuild_rna)
from intervention.operator import (GraphConditionalProcess,response_loss,compatible_pair_loss,
                                   model_contrast_loss)
from intervention.tasks import family,all_actions,control_predictors
from intervention.baselines import SparseResponseSurrogate,KernelResponseSurrogate,break_even
from intervention.metrics import response_error_auc


@pytest.mark.parametrize('task',['edge_switch','node_state'])
def test_valid_exhaustive_family_and_exact_polynomial_controls(task):
    fam=family(task)
    actions=all_actions()
    assert len({fam.transform(a).identity for a in actions})==8
    for order,name in [(1,'additive'),(2,'pairwise'),(3,'higher_order')]:
        f=control_predictors(task)[name]
        base=f(fam.transform())
        y=np.array([f(fam.transform(a))-base for a in actions])
        fit=SparseResponseSurrogate(order=order).condition(actions,y)
        assert np.max(abs(fit.predict(actions)-y))<1e-12
    f=control_predictors(task)['higher_order']
    y=np.array([f(fam.transform(a))-f(fam.transform()) for a in actions])
    fit=SparseResponseSurrogate(order=2).condition(actions,y)
    assert np.max(abs(fit.predict(actions)-y))>.1


def model_episode():
    torch.manual_seed(91226)
    m=GraphConditionalProcess().double()
    fam=family('edge_switch')
    return m,fam,fam.transform(),[fam.transform((0,)),fam.transform((1,))]


def test_node_equivariance_and_graph_invariance():
    m,fam,g,qs=model_episode()
    p=torch.tensor([2,0,3,1])
    def permute(g): return replace(g,x=g.x[p],adjacency=g.adjacency[p][:,p])
    assert torch.allclose(m.edit.graph.node_embeddings(permute(g)),m.edit.graph.node_embeddings(g)[p],atol=1e-12)
    y=torch.tensor([.2,-.3],dtype=torch.float64)
    a=m(g,qs,y,[fam.transform((0,1))])
    b=m(permute(g),[permute(q) for q in qs],y,[permute(fam.transform((0,1)))])
    assert torch.allclose(a,b,atol=1e-12)


def test_transcript_order_identity_and_response_information():
    m,fam,g,qs=model_episode()
    y=torch.tensor([.2,-.3],dtype=torch.float64,requires_grad=True)
    targets=[g,fam.transform((0,1)),fam.transform((0,1,2))]
    p=m(g,qs,y,targets)
    assert p[0].item()==0.
    assert torch.allclose(p,m(g,qs[::-1],y.flip(0),targets),atol=1e-12)
    grad=torch.autograd.grad(p[1:].sum(),y)[0]
    assert grad.abs().sum().item()>1e-10
    assert not torch.allclose(p,m(g,qs,torch.zeros_like(y),targets),atol=1e-10,rtol=0)
    assert not torch.allclose(p,m(g,qs,y.flip(0),targets),atol=1e-10,rtol=0)
    assert m(g,[],torch.empty(0,dtype=torch.float64),[g]).item()==0


def test_no_model_identity_channel_and_indistinguishability():
    m,fam,g,qs=model_episode()
    models=control_predictors('edge_switch')
    f,h=models['indistinguishable_plus'],models['indistinguishable_minus']
    q=[fam.transform(a) for a in all_actions()[:-1]]
    yf=torch.tensor([f(s)-f(g) for s in q],dtype=torch.float64)
    yh=torch.tensor([h(s)-h(g) for s in q],dtype=torch.float64)
    assert torch.equal(yf,yh)
    target=fam.transform((0,1,2))
    a=m(g,q,yf,[target]).item();b=m(g,q,yh,[target]).item()
    assert a==b
    assert max(abs(a-f(target)),abs(b-h(target)))>=abs(f(target)-h(target))/2
    with pytest.raises(TypeError): m(g,q,yf,[target],model_id=42)


def test_auxiliary_loss_backward_does_not_update_weights():
    m,fam,g,qs=model_episode()
    state={k:v.clone() for k,v in m.state_dict().items()}
    y=torch.tensor([.1,.2],dtype=torch.float64)
    targets=[fam.transform(a) for a in [(0,),(1,),(0,1)]]
    p=m(g,qs,y,targets);truth=torch.tensor([.2,-.3,.7],dtype=torch.float64)
    loss=response_loss(p,truth)+compatible_pair_loss(p,truth,[(0,1,2)])
    loss=loss+model_contrast_loss(p,-p,truth,-truth,0.,0.,.01)
    loss.backward()
    assert any(p.grad is not None and p.grad.abs().sum()>0 for p in m.parameters())
    assert all(torch.equal(v,state[k]) for k,v in m.state_dict().items())
    with pytest.raises(ValueError): model_contrast_loss(p,p,truth,truth,0.,1.,.01)


def test_compatibility_rejects_invalid_rectangles():
    fam=family('edge_switch')
    assert fam.combine((0,),(1,))==(0,1)
    with pytest.raises(ValueError): fam.combine((0,),(0,))
    with pytest.raises(ValueError): fam.transform((3,))
    with pytest.raises(ValueError): fam.transform((1,1))


def test_molecule_descriptor_regeneration_identity_and_mapping():
    s='[CH3:1][CH2:2][Cl:3]'
    base=molecular_input(s)
    product,edited,mapping=replace_mapped_atom(s,3,17,9)
    assert mapping=={1:1,2:2,3:3}
    assert base.identity!=edited.identity
    assert not torch.equal(base.descriptors,edited.descriptors)
    assert torch.equal(edited.descriptors,molecular_input(product).descriptors)
    same,identity,_=replace_mapped_atom(s,3,17,17)
    assert same==s and torch.equal(identity.descriptors,base.descriptors)
    assert molecular_input('CCCl').identity==base.identity
    assert molecular_input('F/C=C/F').identity!=molecular_input('F/C=C\\F').identity
    with pytest.raises(ValueError): replace_mapped_atom('CCCl',3,17,9)
    with pytest.raises(ValueError): replace_mapped_atom(s,3,8,9)


def test_molecule_invalid_valence_and_stereo_edit_rejected():
    with pytest.raises(ValueError): replace_mapped_atom('[CH4:1]',1,6,9)
    with pytest.raises(ValueError): replace_mapped_atom('[F:1][C@H:2]([Cl:3])[Br:4]',2,6,7)


def test_registered_library_has_four_corners_or_refuses():
    lib=ProductLibrary({('a','b'):'CC',('x','b'):'CCC',('a','y'):'CO'})
    with pytest.raises(ValueError): lib.combine(('a','b'),{0:'x'},{1:'y'})
    lib=ProductLibrary({('a','b'):'CC',('x','b'):'CCC',('a','y'):'CO',('x','y'):'CCO'})
    assert lib.combine(('a','b'),{0:'x'},{1:'y'})=={0:'x',1:'y'}
    with pytest.raises(ValueError): lib.combine(('a','b'),{0:'x'},{0:'a'})


def test_rna_builder_receives_full_changed_sequence():
    calls=[]
    def builder(s): calls.append(s);return {'sequence':s,'derived':s.count('G')}
    a=rebuild_rna('ACGU',{},builder);b=rebuild_rna('ACGU',{0:'G'},builder)
    assert calls==['ACGU','GCGU'] and b['derived']==a['derived']+1


def test_frozen_queries_cache_every_dependent_input():
    fam=family('edge_switch');g=fam.transform()
    f=lambda s: s.descriptors.sum()+s.context.sum()
    oracle=FrozenQueryOracle(f,'fixture',g.construction)
    assert oracle.effect(g,g)==0 and oracle.query_count==1
    changed=replace(g,descriptors=g.descriptors+1)
    assert oracle.query(changed)!=oracle.query(g) and oracle.query_count==2
    with pytest.raises(ValueError): oracle.query(replace(g,construction='AGILE_Mordred'))
    with pytest.raises(ValueError): oracle.effect(g,replace(g,context=g.context+1))


def test_runtime_auc_counts_cost_and_deadlines():
    truth=np.array([1.,-1.])
    # Exact answer available only after halfway through the log interval.
    assert response_error_auc([.01],[truth],truth,1.,interval=(.001,.1))==pytest.approx(.5)
    assert response_error_auc([.01],[truth],truth,1.,interval=(.001,.1),upfront_s=1000,deployments=1000)==pytest.approx(1.)
    with pytest.raises(ValueError): response_error_auc([.1,.01],[truth,truth],truth,1.)
    assert break_even(10,1,2) is None
    assert break_even(10,2,1)==10


def test_kernel_is_centered_and_conditions_on_queries():
    fit=KernelResponseSurrogate().condition([(),(0,),(1,)],[0.,1.,-.5])
    assert fit.predict([()])[0]==pytest.approx(0.)
    assert np.max(abs(fit.predict([(0,),(1,)])-np.array([1.,-.5])))<1e-5


def test_zero_query_baselines_and_no_completed_output():
    targets=[(),(0,),(0,1,2)]
    for estimator in [SparseResponseSurrogate(),SparseResponseSurrogate(mode='lasso'),KernelResponseSurrogate()]:
        assert np.array_equal(estimator.condition([],[]).predict(targets),np.zeros(3))
        assert estimator.predict([]).shape==(0,)
    assert response_error_auc([],[],[1.,-1.],1.)==pytest.approx(1.)
    with pytest.raises(ValueError): response_error_auc([],[],[1.],1.,deployments=0)
    with pytest.raises(ValueError): response_error_auc([],[],[1.],float('inf'))


def test_known_prior_control_respects_missing_information_and_exact_observations():
    from intervention.baselines import PolynomialGaussianResponseSurrogate,BoostedResponseSurrogate
    gp=PolynomialGaussianResponseSurrogate().condition([(),(0,),(1,),(2,)],[0,1,2,3])
    assert gp.predict([(0,1,2)])[0]==pytest.approx(6.)
    actions=all_actions()
    truth=np.array([0.,1.,2.,3.,5.,8.,10.,20.])
    gp.condition(actions,truth)
    assert np.allclose(gp.predict(actions),truth,atol=1e-12)
    assert np.array_equal(gp.condition([],[]).predict(actions),np.zeros(8))
    assert np.array_equal(BoostedResponseSurrogate().condition([],[]).predict(actions),np.zeros(8))


def test_tabular_features_rebuilt_and_control_factory_does_not_fit():
    from intervention.predictors import fingerprint_descriptor_features,tabular_predictor_control
    base=fingerprint_descriptor_features('CCCl');edited=fingerprint_descriptor_features('CCF')
    assert base.shape==edited.shape==(2052,)
    assert not np.array_equal(base[:2048],edited[:2048])
    assert not np.array_equal(base[-4:],edited[-4:])
    ridge=tabular_predictor_control('ridge')
    assert not hasattr(ridge[-1],'coef_')
    assert not hasattr(tabular_predictor_control('boosted_trees'),'_predictors')
