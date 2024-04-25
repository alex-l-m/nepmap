'''Estimate geometries with XTB'''

from glob import glob
from os.path import basename, splitext, exists
from rdkit.Chem.rdmolfiles import MolFromMolFile, MolToMolFile
from rdkit.Chem.rdmolops import AddHs
from xtb.ase.calculator import XTB
from octahedral_embed import octahedral_embed
from annotate_rdkit_with_ase import optimize_geometry

inpaths = glob('test/*.mol')

for inpath in inpaths:
    mol_id, ext = splitext(basename(inpath))
    outpath = f'xtb_geom/{mol_id}.mol'
    # Skip if it was done already
    if exists(outpath):
        continue
    # Skip the OXA's
    if 'OXA' in mol_id:
        continue
    # Skip the YUY, they're ions
    if 'YUY' in mol_id:
        continue
    # Skip the BIP's, they're ions too
    if 'BIP' in mol_id:
        continue
    mol_nohs = MolFromMolFile(inpath)
    mol = AddHs(mol_nohs)
    mol.SetProp('_Name', mol_id)
    try:
        octahedral_embed(mol, 'tridentate')
    except ValueError:
        print(f'can\'t do {mol_id}')
        continue
    try:
        optimize_geometry(XTB(), mol, conformation_index = 0, uhf = 2)
    except:
        print(f'can\'t do {mol_id}')
        continue
    MolToMolFile(mol, outpath)
