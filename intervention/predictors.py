"""Honest candidate predictor interfaces; no fitting occurs on import."""
import torch
from torch import nn
from .operator import EquivariantGraphEncoder


class FingerprintDescriptorMLP(nn.Module):
    """Tabular MLP control. Fingerprints/descriptors must be freshly regenerated."""
    def __init__(self,input_dim,hidden=64):
        super().__init__()
        self.network=nn.Sequential(nn.Linear(input_dim,hidden),nn.ReLU(),nn.Linear(hidden,1))
    def forward(self,features): return self.network(features).squeeze(-1)


class GraphDescriptorRegressor(nn.Module):
    """Actual message-passing network plus descriptors; not AGILE or TransMA."""
    def __init__(self,node_dim=8,descriptor_dim=4,context_dim=1,hidden=32):
        super().__init__()
        self.encoder=EquivariantGraphEncoder(node_dim,hidden)
        self.readout=nn.Sequential(nn.Linear(hidden+descriptor_dim+context_dim,hidden),
                                   nn.Tanh(),nn.Linear(hidden,1))
    def forward(self,graph):
        return self.readout(torch.cat([self.encoder(graph),graph.descriptors,graph.context])).squeeze(-1)


def fingerprint_descriptor_features(smiles):
    """Morgan radius-2/2048-bit fingerprint plus the four declared RDKit descriptors.

    Regenerate after each edit; no fitted normalization or feature selection.
    This is an explicit candidate representation, not published-model compatibility.
    """
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    from .transforms import molecular_input
    graph=molecular_input(smiles)
    mol=Chem.MolFromSmiles(smiles)
    for atom in mol.GetAtoms():atom.SetAtomMapNum(0)
    fp=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,includeChirality=True).GetFingerprintAsNumPy(mol)
    return np.r_[fp.astype(float),graph.descriptors.numpy()]


def tabular_predictor_control(kind,seed=91226):
    """Unfitted train-only preprocessing and strong predictor control interfaces."""
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    if kind=='ridge':
        from sklearn.linear_model import Ridge
        return make_pipeline(StandardScaler(),Ridge(alpha=1.))
    if kind=='boosted_trees':
        from sklearn.ensemble import HistGradientBoostingRegressor
        return HistGradientBoostingRegressor(max_iter=100,max_leaf_nodes=15,random_state=seed)
    raise ValueError('unknown tabular predictor control')
