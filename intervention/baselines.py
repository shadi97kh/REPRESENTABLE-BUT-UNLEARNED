"""Local response surrogates with a common signed-effect estimand."""
from itertools import combinations
import numpy as np


def binary_design(actions, n=3, order=1):
    terms = [t for k in range(1,order+1) for t in combinations(range(n),k)]
    return np.array([[float(set(t).issubset(a)) for t in terms] for a in actions]).reshape(len(actions),len(terms))


class SparseResponseSurrogate:
    def __init__(self, n=3, order=1, penalty=0., mode='ridge'):
        if mode not in ('ridge','lasso'): raise ValueError('unknown penalty')
        self.n,self.order,self.penalty,self.mode=n,order,penalty,mode

    def condition(self, actions, responses):
        x = binary_design(actions,self.n,self.order)
        y = np.asarray(responses)
        if not len(actions): self.coef=np.zeros(x.shape[1]);return self
        if self.mode == 'lasso':
            from sklearn.linear_model import Lasso
            self.coef = Lasso(alpha=self.penalty,fit_intercept=False,
                              max_iter=10000,tol=1e-8).fit(x,y).coef_
        elif self.penalty == 0:
            self.coef = np.linalg.lstsq(x,y,rcond=None)[0]
        else:
            self.coef = np.linalg.solve(x.T@x+self.penalty*np.eye(x.shape[1]),x.T@y)
        return self

    def predict(self, actions):
        return binary_design(actions,self.n,self.order)@self.coef


class KernelResponseSurrogate:
    """Fixed RBF kernel posterior mean on action indicators, centered at identity."""
    def __init__(self,n=3,lengthscale=1.,noise=1e-6):
        if lengthscale<=0 or noise<=0: raise ValueError('positive kernel parameters required')
        self.n,self.lengthscale,self.noise=n,lengthscale,noise

    def kernel(self,x,y):
        raw = lambda a,b: np.exp(-((a[:,None]-b[None,:])**2).sum(-1)/(2*self.lengthscale**2))
        zero=np.zeros((1,self.n))
        return raw(x,y)-raw(x,zero)-raw(zero,y)+1

    def condition(self,actions,responses):
        self.x=binary_design(actions,self.n,1)
        self.coef=np.linalg.solve(self.kernel(self.x,self.x)+self.noise*np.eye(len(actions)),responses)
        return self

    def predict(self,actions):
        return self.kernel(binary_design(actions,self.n,1),self.x)@self.coef


def break_even(extra_upfront_cost,baseline_online_cost,proposed_online_cost):
    denominator=baseline_online_cost-proposed_online_cost
    return None if denominator<=0 else max(0.,extra_upfront_cost)/denominator


class PolynomialGaussianResponseSurrogate:
    """Noiseless GP mean for the declared iid Gaussian monomial coefficients.

    This is a privileged known-prior control for the analytic fixture family.
    It is not asserted to describe a biological predictor distribution.
    """
    def __init__(self,n=3,order=3): self.n,self.order=n,order
    def condition(self,actions,responses):
        x=binary_design(actions,self.n,self.order)
        self.coefficients=np.linalg.pinv(x)@np.asarray(responses,float)
        return self
    def predict(self,actions):
        return binary_design(actions,self.n,self.order)@self.coefficients


class BoostedResponseSurrogate:
    """Direct local tree response baseline, not a reproduction of ProxySPEX."""
    def __init__(self,n=3,estimators=50,depth=2,learning_rate=.05,seed=91226):
        self.n=n
        self.settings=dict(n_estimators=estimators,max_depth=depth,learning_rate=learning_rate,random_state=seed)
    def condition(self,actions,responses):
        self.model=None
        if actions:
            from sklearn.ensemble import GradientBoostingRegressor
            self.model=GradientBoostingRegressor(**self.settings).fit(binary_design(actions,self.n),responses)
        return self
    def predict(self,actions):
        if not actions:return np.empty(0)
        if self.model is None:return np.zeros(len(actions))
        return self.model.predict(binary_design(actions,self.n))-self.model.predict(np.zeros((1,self.n)))[0]
