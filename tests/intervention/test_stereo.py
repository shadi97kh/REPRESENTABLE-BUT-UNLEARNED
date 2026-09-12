import torch
from rdkit import Chem
from intervention.transforms import molecular_input
from intervention.operator import EquivariantGraphEncoder


def test_stereo_encoder_preserves_isomer_and_smiles_order():
    torch.manual_seed(91226)
    encoder=EquivariantGraphEncoder().double()
    molecule=Chem.MolFromSmiles('C[C@H](F)CCl')
    a=molecular_input(Chem.MolToSmiles(molecule,canonical=False,rootedAtAtom=0))
    b=molecular_input(Chem.MolToSmiles(molecule,canonical=False,rootedAtAtom=4))
    assert a.identity==b.identity
    assert torch.allclose(encoder(a),encoder(b),atol=1e-12)
    e=molecular_input('F/C=C/Cl');z=molecular_input('F/C=C\\Cl')
    assert not torch.equal(e.edge_stereo,z.edge_stereo)
    assert not torch.allclose(encoder(e),encoder(z),atol=1e-12)
