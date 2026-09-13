"""Conventional compact learned message passing and matched controls; no novelty assertion."""
import os
os.environ.setdefault('CUBLAS_WORKSPACE_CONFIG',':4096:8')
import numpy as np
import torch
from torch import nn

class GraphModel(nn.Module):
 def __init__(self,features,context,nodes=64,hidden=32,layers=2,kind='gnn',dropout=.1,chem_start=39,chem_end=100):
  super().__init__();self.kind=kind;self.chem_start=chem_start;self.chem_end=chem_end;self.hidden=hidden
  self.embed=nn.Linear(features,hidden);self.self_layers=nn.ModuleList([nn.Linear(hidden,hidden) for _ in range(layers)]);self.norms=nn.ModuleList([nn.LayerNorm(hidden) for _ in range(layers)])
  self.message_weights=nn.ParameterList([nn.Parameter(torch.empty(8,hidden,hidden)) for _ in range(layers)])
  for w in self.message_weights:nn.init.xavier_uniform_(w)
  self.dropout=nn.Dropout(dropout);self.readout=nn.Sequential(nn.Linear(nodes*hidden+context,64),nn.ReLU(),nn.Dropout(dropout),nn.Linear(64,1))
 def forward(self,X,A,mask,C):
  if self.kind=='gnn_no_chemistry':
   X=X.clone();X[:,:,self.chem_start:self.chem_end]=0
   A=A.clone();A[:,0]+=A[:,1]+A[:,2];A[:,3]+=A[:,4]+A[:,5];A[:,[1,2,4,5]]=0
  h=torch.relu(self.embed(X))*mask[:,:,None];degree=A.sum((1,3)).clamp(min=1)[:,:,None]
  for own,w,norm in zip(self.self_layers,self.message_weights,self.norms):
   if self.kind=='nongraph':agg=A.sum(-1)[:,:,:,None]*h[:,None,:,:]
   else:agg=torch.matmul(A,h[:,None,:,:])
   message=torch.einsum('brnh,rhk->bnk',agg,w)/degree
   h=norm(h+self.dropout(torch.relu(own(h)+message)))*mask[:,:,None]
  return self.readout(torch.cat([h.flatten(1),C],1)).squeeze(-1)

class TokenCNN(nn.Module):
 def __init__(self,features,context,nodes=64,hidden=32,**kwargs):
  super().__init__();self.embed=nn.Linear(features,hidden);self.layers=nn.Sequential(nn.Conv1d(hidden,hidden,3,padding=1),nn.ReLU(),nn.Dropout(.1),nn.Conv1d(hidden,hidden,3,padding=1),nn.ReLU());self.readout=nn.Sequential(nn.Linear(nodes*hidden+context,64),nn.ReLU(),nn.Dropout(.1),nn.Linear(64,1))
 def forward(self,X,A,mask,C):
  h=torch.relu(self.embed(X))*mask[:,:,None];h=(h+self.layers(h.transpose(1,2)).transpose(1,2))*mask[:,:,None]
  return self.readout(torch.cat([h.flatten(1),C],1)).squeeze(-1)

def build(kind,manifest,config):
 kwargs=dict(features=manifest['node_features'],context=len(manifest['context_features']),nodes=manifest['max_nodes'],hidden=config['hidden'])
 if kind=='token_cnn':return TokenCNN(**kwargs)
 # No-chemistry removes modification, stereochemistry/terminal metadata but preserves position/end markers by an explicit input mask below.
 return GraphModel(**kwargs,kind=kind,chem_start=manifest['chemistry_columns'][0],chem_end=manifest['node_features'])

def tensors(path,device):
 z=np.load(path);return {k:torch.as_tensor(z[k],device=device) for k in ['X','A','mask','C']}

def take(data,idx):return {k:v[idx] for k,v in data.items()}

def predict(model,data,batch=128):
 model.eval();out=[]
 with torch.no_grad():
  for i in range(0,len(data['X']),batch):out.append(model(**take(data,slice(i,i+batch))).cpu().numpy())
 return np.concatenate(out)
