"""Predeclared primary metric and paired independent-function uncertainty."""
import numpy as np


def response_error_auc(completed_times, completed_predictions, truth, scale,
                       interval=(.001,1.), upfront_s=0., deployments=1000):
    """Exact integral of normalized MSE of the latest completed output vs log time.

    Zero-effect predictions apply before the first completed output; all actions
    receive a prediction. Late/timeout runs retain their last completed output.
    No optimistic interpolation or exclusion of hard actions is permitted.
    """
    lo,hi=interval
    if not (0<lo<hi and np.isfinite(scale) and scale>0 and deployments>0 and np.isfinite(upfront_s) and upfront_s>=0): raise ValueError('invalid metric parameters')
    times=np.asarray(completed_times,float)+upfront_s/deployments
    preds=np.asarray(completed_predictions,float)
    truth=np.asarray(truth,float)
    if truth.ndim!=1 or not len(truth): raise ValueError('nonempty one-dimensional target required')
    if not len(times) and not preds.size: preds=np.empty((0,len(truth)))
    if times.ndim!=1 or len(times)!=len(preds) or (np.diff(times)<0).any() or (times<0).any():
        raise ValueError('times must be nonnegative, ordered completions')
    if preds.shape != (len(times),len(truth)): raise ValueError('predictions must cover every action')
    if not all(np.isfinite(x).all() for x in (times,preds,truth)): raise ValueError('nonfinite curve')
    knots=np.unique(np.r_[lo,times[(times>lo)&(times<hi)],hi])
    area=0.
    for left,right in zip(knots[:-1],knots[1:]):
        idx=np.searchsorted(times,left,side='right')-1
        p=np.zeros_like(truth) if idx<0 else preds[idx]
        area+=np.mean(((p-truth)/scale)**2)*np.log(right/left)
    return float(area/np.log(hi/lo))


def paired_function_interval(proposed,baseline,seed=91226,replicates=10000):
    """Inputs are per-independent-function AUCs, already averaged across graphs."""
    p,b=np.asarray(proposed,float),np.asarray(baseline,float)
    if p.ndim!=1 or p.shape!=b.shape or len(p)<2: raise ValueError('matched independent function units required')
    if not np.isfinite(p).all() or not np.isfinite(b).all(): raise ValueError('nonfinite results')
    d=b-p
    rng=np.random.default_rng(seed)
    means=d[rng.integers(len(d),size=(replicates,len(d)))].mean(1)
    return dict(improvement=float(d.mean()),ci95=np.quantile(means,[.025,.975]).tolist(),
                relative_improvement=None if b.mean()<=0 else float(d.mean()/b.mean()),n_functions=len(d))
