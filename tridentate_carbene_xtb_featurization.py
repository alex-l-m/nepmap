from glob import glob
from os.path import basename, splitext
import csv
from rdkit.Chem.rdmolfiles import MolFromSmarts, MolFromMolFile
from rdkit.Chem.Draw import MolToImageFile
from ase import Atoms
from ase.optimize import BFGS
from ase.units import Bohr
from xtb.ase.calculator import XTB
from xtb.interface import Calculator, Param
from annotate_rdkit_with_ase import rdkit2ase

pattern = MolFromSmarts('[#7]~[#6](~[Ir]12)~[#7]~[#6]~[#6]~1~[#6]~[#7]~[#6]~2~[#7]')
MolToImageFile(pattern, 'featurization_pattern.png')

tridentate_carbene_paths = glob('xtb_geom/*.mol')

output_featurization_path = 'featurization.csv'
# Initialize file with header
header = ['mol_id', 'feature_id', 'feature_value']
with open(output_featurization_path, 'w') as f:
    writer = csv.writer(f)
    writer.writerow(header)

for path in tridentate_carbene_paths:
    # Load tridentate carbene as an rdkit molecule
    mol_id, _ = splitext(basename(path))
    rdkit_mol = MolFromMolFile(path, removeHs = False)
    if rdkit_mol is None:
        continue
    has_pattern = rdkit_mol.HasSubstructMatch(pattern)
    assert has_pattern

    target_atom_indices = rdkit_mol.GetSubstructMatch(pattern)

    target_bond_index_pairs = []
    for bond in pattern.GetBonds():
        begin_idx = target_atom_indices[bond.GetBeginAtomIdx()]
        end_idx = target_atom_indices[bond.GetEndAtomIdx()]
        target_bond_index_pairs.append((begin_idx, end_idx))

    # Sanity check: all bonds connect atoms in the pattern, and are actually
    # bonded in the molecule
    for begin_index, end_index in target_bond_index_pairs:
        assert begin_index in target_atom_indices
        assert end_index in target_atom_indices
        assert rdkit_mol.GetBondBetweenAtoms(begin_index, end_index) is not None

    # Convert to ASE atoms
    ase_atoms = rdkit2ase(rdkit_mol, conformation_index=0, charge=0, uhf=2)
    ase_atoms.set_calculator(XTB(method='GFN2-xTB'))

    # Do XTB single-point calculation, also in triplet state
    # Lots of copying from docs example:
    # https://xtb-python.readthedocs.io/en/latest/general-api.html#xtb.interface.Calculator
    # Imitating how XTB ASE calculator does the conversion:
    # https://github.com/grimme-lab/xtb-python/blob/main/xtb/ase/calculator.py
    positions = ase_atoms.get_positions() / Bohr
    numbers = ase_atoms.get_atomic_numbers()
    calc = Calculator(Param.GFN2xTB, numbers, positions, uhf=2)
    res = calc.singlepoint()

    # Extract charges and bond orders of pattern from XTB results
    all_charges = res.get_charges()
    pattern_charges = [all_charges[i] for i in target_atom_indices]
    all_bond_orders = res.get_bond_orders()
    pattern_bond_orders = [all_bond_orders[i, j] \
            for i, j in target_bond_index_pairs]

    # Write results by appending to output featurization file
    with open(output_featurization_path, 'a') as f:
        writer = csv.writer(f)
        for i, charge in enumerate(pattern_charges):
            writer.writerow([mol_id, f'charge_{i}', charge])
        for i, bond_order in enumerate(pattern_bond_orders):
            writer.writerow([mol_id, f'bond_order_{i}', bond_order])
