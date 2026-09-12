"""Complete input reconstruction for valid interventions, independent of scmp."""
from __future__ import annotations
from dataclasses import dataclass, replace
import hashlib
import json
from typing import Callable
import numpy as np
import torch


@dataclass(frozen=True)
class GraphInput:
    x: torch.Tensor
    adjacency: torch.Tensor
    descriptors: torch.Tensor
    context: torch.Tensor
    construction: str
    identity: str
    edge_stereo: torch.Tensor | None = None

    def validate(self):
        n = len(self.x)
        if self.x.ndim != 2 or not n or self.adjacency.shape != (n, n):
            raise ValueError('invalid graph dimensions')
        if self.descriptors.ndim != 1 or self.context.ndim != 1:
            raise ValueError('descriptors/context must be vectors')
        if not all(torch.isfinite(x).all() for x in
                   (self.x, self.adjacency, self.descriptors, self.context)):
            raise ValueError('nonfinite predictor input')
        if not torch.equal(self.adjacency, self.adjacency.T) or (self.adjacency < 0).any():
            raise ValueError('expected undirected nonnegative adjacency')
        if torch.any(self.adjacency.diag() != 0):
            raise ValueError('self loops are an encoder operation, not input edges')
        if self.edge_stereo is not None:
            if self.edge_stereo.shape != (n,n,2) or not torch.isfinite(self.edge_stereo).all():
                raise ValueError("invalid edge stereochemistry")
            if not torch.equal(self.edge_stereo,self.edge_stereo.transpose(0,1)):
                raise ValueError("edge stereochemistry must be symmetric")
        return self


def digest(parts):
    return hashlib.sha256(json.dumps(parts, sort_keys=True).encode()).hexdigest()


@dataclass(frozen=True)
class BinaryGraphFamily:
    """Three declared, commuting switches; all eight states are valid."""
    task: str
    x: torch.Tensor
    adjacency: torch.Tensor
    switches: tuple

    def transform(self, action=()):
        if len(set(action)) != len(action) or any(i not in range(len(self.switches)) for i in action):
            raise ValueError('actions are sets of distinct declared switches')
        x, adj = self.x.clone(), self.adjacency.clone()
        for i in action:
            if self.task == 'edge_switch':
                a,b = self.switches[i]
                adj[a,b] = adj[b,a] = 1 - adj[a,b]
            elif self.task == 'node_state':
                node = self.switches[i]
                x[node,0] = 1 - x[node,0]
            else:
                raise ValueError('undeclared task')
        # All derived features are reconstructed after every change.
        desc = torch.stack([adj.sum()/2, x[:,0].sum(), (adj@adj@adj).trace()/6,
                            (adj.sum(1)**2).sum()]).to(x)
        ident = digest(dict(task=self.task, x=x.tolist(), adjacency=adj.tolist()))
        return GraphInput(x,adj,desc,torch.zeros(1,dtype=x.dtype),
                          'tiny_'+self.task+'_v1',ident).validate()

    def combine(self, a, b):
        if set(a) & set(b):
            raise ValueError('overlapping switches are not a compatible interaction pair')
        result = tuple(sorted((*a,*b)))
        for state in [(),a,b,result]: self.transform(state)
        return result


MOLECULAR_SCHEMA = 'rdkit_graph8_desc4_stereo2_v2; NOT AGILE Mordred or TransMA 3D'


def molecular_input(smiles, context=(0.,)):
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('invalid molecule')
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    x = [[a.GetAtomicNum()/100, a.GetFormalCharge(), a.GetIsAromatic(),
          a.GetTotalNumHs(), a.GetDegree(), ({"R":1,"S":2}.get(a.GetProp("_CIPCode"),0) if a.HasProp("_CIPCode") else 0),
          a.GetMass()/100, a.IsInRing()] for a in mol.GetAtoms()]
    adj = np.zeros((len(x),len(x)))
    stereo = np.zeros((len(x),len(x),2))
    for b in mol.GetBonds():
        i,j = b.GetBeginAtomIdx(),b.GetEndAtomIdx()
        adj[i,j] = adj[j,i] = b.GetBondTypeAsDouble()
        for channel,tag in enumerate((Chem.BondStereo.STEREOE,Chem.BondStereo.STEREOZ)):
            stereo[i,j,channel] = stereo[j,i,channel] = float(b.GetStereo()==tag)
    desc = [Descriptors.MolWt(mol), Crippen.MolLogP(mol),
            rdMolDescriptors.CalcTPSA(mol), rdMolDescriptors.CalcNumRotatableBonds(mol)]
    # Atom-map labels are correspondence metadata, not molecular identity.
    keymol = Chem.Mol(mol)
    for a in keymol.GetAtoms(): a.SetAtomMapNum(0)
    identity = Chem.MolToSmiles(keymol, isomericSmiles=True, canonical=True)
    t = lambda a: torch.as_tensor(a,dtype=torch.float64)
    return GraphInput(t(x),t(adj),t(desc),t(context),MOLECULAR_SCHEMA,identity,t(stereo)).validate()


def replace_mapped_atom(smiles, map_id, expected_atomic_num, new_atomic_num, context=(0.,)):
    """Small supported edit: mapped atom substitution, fixed connectivity.

    Reject edits at specified tetrahedral stereocenters. Stable atom maps record
    correspondence; no canonical-string position is treated as an atom ID.
    More general fragment transformations require an explicit registered product.
    """
    from rdkit import Chem
    mol = Chem.MolFromSmiles(smiles)
    if mol is None: raise ValueError('invalid reference')
    ids = [a.GetAtomMapNum() for a in mol.GetAtoms()]
    if any(i <= 0 for i in ids) or len(set(ids)) != len(ids):
        raise ValueError('complete unique atom maps required')
    if map_id not in ids: raise ValueError('unknown atom map')
    atom = mol.GetAtomWithIdx(ids.index(map_id))
    if atom.GetAtomicNum() != expected_atomic_num: raise ValueError('reference atom mismatch')
    if new_atomic_num == expected_atomic_num:
        return smiles, molecular_input(smiles, context), {i:i for i in ids}
    if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED:
        raise ValueError('stereocenter substitution requires a declared stereochemical product')
    atom.SetAtomicNum(new_atomic_num)
    try: Chem.SanitizeMol(mol)
    except Exception as e: raise ValueError('invalid edited valence') from e
    product = Chem.MolToSmiles(mol,isomericSmiles=True)
    return product, molecular_input(product, context), {i:i for i in ids}


class ProductLibrary:
    """Component substitutions resolve to observed, explicitly registered SMILES.

    Chemical validity does not establish assay comparability or synthesis yield.
    This interface makes no atom correspondence claim for fragment substitutions.
    """
    def __init__(self, products: dict[tuple, str]):
        self.products = dict(products)

    def transform(self, reference, substitutions=None, context=(0.,)):
        key = list(reference)
        for axis,value in (substitutions or {}).items():
            if not 0 <= axis < len(key): raise ValueError('invalid component axis')
            key[axis] = value
        key = tuple(key)
        if key not in self.products: raise ValueError('unregistered product')
        return molecular_input(self.products[key],context)

    def combine(self, reference, a, b):
        if a.keys() & b.keys(): raise ValueError('overlapping component edits')
        ab = {**a, **b}
        for s in ({},a,b,ab): self.transform(reference,s)
        return ab


def rebuild_rna(sequence, edits: dict[int,str], feature_builder: Callable):
    """The registered builder receives the complete edited RNA sequence.

    It must regenerate every sequence/structure feature and expose its version.
    No isolated structural or biological-intervention claim is made.
    """
    if not sequence or set(sequence)-set('ACGU'): raise ValueError('invalid RNA')
    out = list(sequence)
    for i,b in edits.items():
        if i not in range(len(out)) or b not in ('A','C','G','U'): raise ValueError('invalid RNA edit')
        out[i] = b
    return feature_builder(''.join(out))


class FrozenQueryOracle:
    """Counts uncached f evaluations; cache keys include construction and context."""
    def __init__(self, predictor, model_hash, required_construction):
        self.predictor, self.model_hash = predictor, model_hash
        self.required_construction = required_construction
        self.cache, self.query_count = {}, 0
        if isinstance(predictor, torch.nn.Module):
            predictor.eval()
            for p in predictor.parameters(): p.requires_grad_(False)

    def query(self, state):
        state.validate()
        if state.construction != self.required_construction:
            raise ValueError('predictor feature construction mismatch')
        key = digest(dict(model=self.model_hash,construction=state.construction,
                          identity=state.identity,context=state.context.tolist(),
                          x=state.x.tolist(),adjacency=state.adjacency.tolist(),
                          descriptors=state.descriptors.tolist(),
                          edge_stereo=None if state.edge_stereo is None else state.edge_stereo.tolist()))
        if key not in self.cache:
            with torch.no_grad(): value = float(self.predictor(state))
            if not np.isfinite(value): raise ValueError('nonfinite query response')
            self.cache[key] = value
            self.query_count += 1
        return self.cache[key]

    def effect(self, reference, edited):
        if not torch.equal(reference.context,edited.context):
            raise ValueError('context changed in a fixed-context intervention')
        return self.query(edited)-self.query(reference)


def interaction(oracle, reference, a, b, ab):
    """Caller must obtain four states through a compatible-action constructor."""
    return oracle.effect(reference,ab)-oracle.effect(reference,a)-oracle.effect(reference,b)
