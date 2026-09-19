"""An identity-centered graph conditional neural process, not a novelty claim."""
from __future__ import annotations
import torch
from torch import nn


class EquivariantGraphEncoder(nn.Module):
    def __init__(self, node_dim=8, hidden=32, layers=2):
        super().__init__()
        self.input = nn.Linear(node_dim,hidden)
        self.self_layers = nn.ModuleList(nn.Linear(hidden,hidden) for _ in range(layers))
        self.neighbor_layers = nn.ModuleList(nn.Linear(hidden,hidden,bias=False) for _ in range(layers))
        self.stereo_layers = nn.ModuleList(nn.ModuleList(nn.Linear(hidden,hidden,bias=False) for _ in range(2)) for _ in range(layers))

    def node_embeddings(self, graph):
        h = torch.tanh(self.input(graph.x))
        # Sum aggregation preserves node-permutation equivariance.
        for own,neighbor,stereo in zip(self.self_layers,self.neighbor_layers,self.stereo_layers):
            messages = neighbor(graph.adjacency@h)
            if graph.edge_stereo is not None:
                messages = messages + sum(layer(graph.edge_stereo[:,:,j]@h) for j,layer in enumerate(stereo))
            h = torch.tanh(own(h)+messages)
        return h

    def forward(self, graph):
        return self.node_embeddings(graph).sum(0)


class GraphEditEncoder(nn.Module):
    def __init__(self, node_dim=8, descriptor_dim=4, context_dim=1, hidden=32):
        super().__init__()
        self.graph = EquivariantGraphEncoder(node_dim,hidden)
        self.projection = nn.Sequential(nn.Linear(hidden*3+descriptor_dim*3+context_dim,hidden),
                                        nn.Tanh(),nn.Linear(hidden,hidden))

    def forward(self, base, edited):
        if not torch.equal(base.context,edited.context): raise ValueError('fixed context required')
        b,a = self.graph(base),self.graph(edited)
        x = torch.cat([b,a,a-b,base.descriptors,edited.descriptors,
                       edited.descriptors-base.descriptors,base.context])
        return self.projection(x)


class GraphConditionalProcess(nn.Module):
    """Inputs contain actions and observed responses; no model ID/assay labels.

    The response-only proposed operator and capacity-matched generic CNP are the
    SAME class. Auxiliary losses are separate ablations, not a new architecture.
    """
    def __init__(self, hidden=32, node_dim=8, descriptor_dim=4, context_dim=1):
        super().__init__()
        self.edit = GraphEditEncoder(node_dim,descriptor_dim,context_dim,hidden)
        self.set_encoder = nn.Sequential(nn.Linear(hidden+1,hidden),nn.Tanh(),nn.Linear(hidden,hidden))
        self.decoder = nn.Sequential(nn.Linear(2*hidden+1,hidden),nn.Tanh(),nn.Linear(hidden,1))
        self.hidden = hidden

    def encode_transcript(self, base, query_states, responses):
        if len(query_states) != len(responses): raise ValueError('query/response length mismatch')
        if not torch.isfinite(responses).all(): raise ValueError('nonfinite transcript')
        if not query_states:
            return base.x.new_zeros(self.hidden),base.x.new_zeros(1)
        encoded = torch.stack([self.edit(base,g) for g in query_states])
        z = self.set_encoder(torch.cat([encoded,responses.reshape(-1,1)],dim=-1)).mean(0)
        count = base.x.new_tensor([len(query_states)]).log1p()
        return z,count

    def forward(self, base, query_states, responses, targets):
        z,count = self.encode_transcript(base,query_states,responses)
        def decode(g): return self.decoder(torch.cat([self.edit(base,g),z,count])).squeeze(-1)
        identity = decode(base)
        if not targets: return base.x.new_empty(0)
        return torch.stack([decode(g)-identity for g in targets])


def response_loss(predictions, responses, scale=1.):
    if scale <= 0: raise ValueError('positive train-derived scale required')
    return ((predictions-responses)/scale).square().mean()


def compatible_pair_loss(predictions, responses, triples, scale=1.):
    """Triples index compatible (a,b,ab) states already checked by a family."""
    if not triples: return predictions.sum()*0
    idx = torch.tensor(triples,dtype=torch.long,device=predictions.device)
    def second(x): return x[idx[:,2]]-x[idx[:,0]]-x[idx[:,1]]
    return response_loss(second(predictions),second(responses),scale)


def model_contrast_loss(pred_f,pred_g,response_f,response_g,base_f,base_g,tolerance,scale=1.):
    """Source/development model pairs only; no test model-selection path."""
    if abs(base_f-base_g)>tolerance: raise ValueError('not a prespecified matched-output pair')
    return response_loss(pred_f-pred_g,response_f-response_g,scale)
