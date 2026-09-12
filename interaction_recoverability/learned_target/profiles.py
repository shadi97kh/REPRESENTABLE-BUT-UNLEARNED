"""Globally constrained analytic profiles, float64 CPU only."""
import math
import numpy as np
import torch
from torch import nn

torch.set_num_threads(1)
try:torch.set_num_interop_threads(1)
except RuntimeError:pass
torch.set_default_dtype(torch.float64)

ACTIONS=np.vstack([np.zeros(4),np.eye(4),[1,1,0,0]])
FEATURES=[('sin',w,0.) for w in [.5,1.,1.5,2.,2.5,3.]]+[
    ('gauss',.75,c) for c in [-3.,-2.,-1.,0.,1.,2.,3.]]+[
    ('tanh',1.,c) for c in [-2.,0.,2.]]
K=len(FEATURES)

def tensor(x):return torch.as_tensor(x,dtype=torch.float64)

def warp(z,delta):
    return z+delta*z*torch.sqrt(1+z*z)

def warp_derivatives(z,delta):
    hp=1+delta*(1+2*z*z)/torch.sqrt(1+z*z)
    hpp=delta*z*(3+2*z*z)/(1+z*z)**1.5
    return hp,hpp

def atoms(x,max_order=3):
    """Vectorized analytic atoms; global derivative bounds (2,3,9,27)."""
    xx=x[...,None]
    frequencies=tensor([.5,1.,1.5,2.,2.5,3.])
    centers=tensor([-3.,-2.,-1.,0.,1.,2.,3.])
    u=.75*(xx-centers);eg=torch.exp(-u*u)
    t=torch.tanh(xx-tensor([-2.,0.,2.]));d1=1-t*t
    gaussian=[eg-torch.exp(-(.75*centers)**2)]
    hyper=[t-torch.tanh(-tensor([-2.,0.,2.]))]
    if max_order>=1:
        gaussian.append(-1.5*u*eg);hyper.append(d1)
    if max_order>=2:
        gaussian.append(.75**2*(4*u*u-2)*eg);hyper.append(-2*t*d1)
    if max_order>=3:
        gaussian.append(.75**3*(-8*u**3+12*u)*eg);hyper.append(d1*(-2+6*t*t))
    return [torch.cat([frequencies**d*torch.sin(xx*frequencies+d*math.pi/2),gaussian[d],hyper[d]],-1) for d in range(max_order+1)]


def relative_profile_directions(x):
    k=atoms(x)
    return [k[0],k[0]+k[1],k[0]+2*k[1]+k[2],k[0]+3*k[1]+3*k[2]+k[3]]

def project_l1(x,radius=1.):
    if float(x.abs().sum())<=radius:return x
    u=torch.sort(x.abs(),descending=True).values
    css=torch.cumsum(u,0)-radius
    index=torch.arange(1,len(u)+1,dtype=x.dtype)
    rho=torch.nonzero(u-css/index>0)[-1,0]
    theta=css[rho]/(rho+1)
    return x.sign()*torch.clamp(x.abs()-theta,min=0)

class ProfileModel(nn.Module):
    def __init__(self,kind='neural',width=4,seed=0):
        super().__init__();self.kind=kind;self.width=width
        rng=torch.Generator().manual_seed(int(seed))
        self.coef=nn.Parameter(.75+.01*torch.randn((2,2),generator=rng))
        if kind=='neural':
            self.a=nn.Parameter(.05*torch.randn((2,width),generator=rng))
            self.w=nn.Parameter(.25+.75*torch.rand((2,width),generator=rng))
            self.b=nn.Parameter(-1.5+3*torch.rand((2,width),generator=rng))
        elif kind=='basis':self.theta=nn.Parameter(torch.zeros((2,K)))
        else:raise ValueError('neural or basis required')
        self.project()

    @torch.no_grad()
    def project(self):
        self.coef.clamp_(.501,.999)
        if self.kind=='neural':
            self.a.clamp_(-1,1);self.w.clamp_(.25,1);self.b.clamp_(-3,3)
        else:
            for j in range(2):self.theta[j].copy_(project_l1(self.theta[j]))

    def profile(self,j,x,order=0):
        if self.kind=='neural':
            t=torch.tanh(x[...,None]*self.w[j]+self.b[j]);d1=1-t*t
            terms=[t-torch.tanh(self.b[j]),self.w[j]*d1,
                   -2*self.w[j]**2*t*d1,self.w[j]**3*d1*(-2+6*t*t)]
            mod=sum(math.comb(order,k)*(terms[k]*self.a[j]).mean(-1) for k in range(order+1))
            return torch.exp(x)*(1+.1*mod)
        features=atoms(x,order)
        rel=sum(math.comb(order,k)*features[k] for k in range(order+1))
        return torch.exp(x)*(1+.08*(rel@self.theta[j]))

    def response(self,a,z,delta,derivatives=False):
        aa=tensor(ACTIONS[a] if isinstance(a,int) else a)
        x=z+self.coef[0,0]*aa[0]+self.coef[1,0]*aa[1]+aa[2]
        y=warp(z,delta)+self.coef[0,1]*aa[0]+self.coef[1,1]*aa[1]+aa[3]
        F=self.profile(0,x)+self.profile(1,y)
        if not derivatives:return F
        hp,hpp=warp_derivatives(z,delta)
        Fp=self.profile(0,x,1)+self.profile(1,y,1)*hp
        Fpp=self.profile(0,x,2)+self.profile(1,y,2)*hp*hp+self.profile(1,y,1)*hpp
        return F,Fp,Fpp,x,y,hp

    def payload(self):
        return {'kind':self.kind,'width':self.width,'state':{k:v.detach().tolist() for k,v in self.state_dict().items()}}

    @classmethod
    def from_payload(cls,payload):
        model=cls(payload['kind'],payload['width'])
        model.load_state_dict({k:tensor(v) for k,v in payload['state'].items()})
        return model

def inverse(model,a,y,delta,implicit=False):
    """Bracketed inverse; optional exact first implicit derivative at fixed outcome."""
    y=tensor(y)
    if not torch.isfinite(y).all() or not (y>0).all():raise ValueError('positive finite outcomes required')
    with torch.no_grad():
        lo=torch.full_like(y,-24.);hi=torch.full_like(y,24.)
        if (model.response(a,lo,delta)>y).any() or (model.response(a,hi,delta)<y).any():
            raise ArithmeticError('outcome outside declared inverse bracket [-24,24]')
        for _ in range(52):
            mid=(lo+hi)/2;below=model.response(a,mid,delta)<y
            lo=torch.where(below,mid,lo);hi=torch.where(below,hi,mid)
        z=(lo+hi)/2
        residual=(model.response(a,z,delta)-y).abs()/(1+y)
        if float(residual.max())>1e-10:raise ArithmeticError('inverse residual exceeds tolerance')
    if implicit:
        F,Fp,*_=model.response(a,z,delta,True)
        z=z+(y-F)/Fp.detach()
    return z

def negative_log_likelihood(model,a,y,delta):
    z=inverse(model,a,y,delta,implicit=True)
    _,Fp,*_=model.response(a,z,delta,True)
    return .5*z*z+torch.log(Fp)+.5*math.log(2*math.pi)

def quadrature(nodes=96,radius=10.):
    x,w=np.polynomial.legendre.leggauss(nodes);z=radius*x
    weights=radius*w*np.exp(-z*z/2)/math.sqrt(2*math.pi)
    weights/=weights.sum()
    return tensor(z),tensor(weights)

def target(model,delta,nodes=96):
    z,w=quadrature(nodes)
    return w@(model.response(5,z,delta)-model.response(1,z,delta)-model.response(2,z,delta)+model.response(0,z,delta))
