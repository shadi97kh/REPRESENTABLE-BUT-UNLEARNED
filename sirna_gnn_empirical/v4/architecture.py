"""Unchanged v2 conventional architecture. Correction affects relation sharing only.
No-message uses each receiver state times its relation degree; it has no neighbor-state dependence.
"""
import torch
from torch import nn
class GraphModel(nn.Module):
 def __init__(self,kind='gnn',R=False,support=None):
  super().__init__();self.kind=kind;self.R=R;self.embed=nn.Linear(93,32);self.self_layers=nn.ModuleList([nn.Linear(32,32) for _ in range(2)]);self.norms=nn.ModuleList([nn.LayerNorm(32) for _ in range(2)]);self.message_weights=nn.ParameterList([nn.Parameter(torch.empty(8,32,32)) for _ in range(2)])
  for w in self.message_weights:nn.init.xavier_uniform_(w)
  self.dropout=nn.Dropout(.1);self.readout=nn.Sequential(nn.Linear(64*32+8,64),nn.ReLU(),nn.Dropout(.1),nn.Linear(64,1));self.support=support or [True]*8
 def forward(self,X,A,mask,C):
  if self.kind=='gnn_no_chemistry':
   X=X.clone();X[:,:,39:87]=0;X[:,:,89:93]=0
   A=A.clone();A[:,2]+=A[:,0]+A[:,1];A[:,5]+=A[:,3]+A[:,4];A[:,[0,1,3,4]]=0
  h=torch.relu(self.embed(X))*mask[:,:,None];degree=A.sum((1,3)).clamp(min=1)[:,:,None]
  for own,w,norm in zip(self.self_layers,self.message_weights,self.norms):
   if self.R:
    w=torch.stack([w[2]+w[0]*self.support[0],w[2]+w[1]*self.support[1],w[2],w[5]+w[3]*self.support[3],w[5]+w[4]*self.support[4],w[5],w[6]*self.support[6],w[7]*self.support[7]])
   agg=A.sum(-1)[:,:,:,None]*h[:,None,:,:] if self.kind=='nongraph' else torch.matmul(A,h[:,None,:,:])
   message=torch.einsum('brnh,rhk->bnk',agg,w)/degree
   h=norm(h+self.dropout(torch.relu(own(h)+message)))*mask[:,:,None]
  return self.readout(torch.cat([h.flatten(1),C],1)).squeeze(-1)
class TokenCNN(nn.Module):
 def __init__(self):
  super().__init__();self.embed=nn.Linear(93,32);self.layers=nn.Sequential(nn.Conv1d(32,32,3,padding=1),nn.ReLU(),nn.Dropout(.1),nn.Conv1d(32,32,3,padding=1),nn.ReLU());self.readout=nn.Sequential(nn.Linear(2056,64),nn.ReLU(),nn.Dropout(.1),nn.Linear(64,1))
 def forward(self,X,A,mask,C):
  h=torch.relu(self.embed(X))*mask[:,:,None];h=(h+self.layers(h.transpose(1,2)).transpose(1,2))*mask[:,:,None];return self.readout(torch.cat([h.flatten(1),C],1)).squeeze(-1)
