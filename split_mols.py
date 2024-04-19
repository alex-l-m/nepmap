'''Load molecules from mol2 files, split into individual ligands, and output smiles strings of ligands to a csv file'''
from glob import glob
from os.path import basename, splitext
import pandas as pd
from rdkit.Chem.rdmolfiles import MolFromMol2File
from rdkit.Chem.rdmolops import RemoveStereochemistry
from modify_organometallics import adjust_iridium_bonds, fragment

inpaths = glob('mol2/*.mol2')

rows = []

for chemdraw in [True, False]:
    for inpath in inpaths:
        mol_id, ext = splitext(basename(inpath))
        raw_mol = MolFromMol2File(inpath)
        # ELOYIP doesn't sanitize
        if raw_mol is None:
            ligand_a = None
            ligand_b = None
        else:
            RemoveStereochemistry(raw_mol)
            mol = adjust_iridium_bonds(raw_mol, chemdraw = chemdraw)
            # Set formal charge on iridium to 0
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "Ir":
                    atom.SetFormalCharge(0)
            ligands = fragment(mol)
            # Keep two largest ones
            # Needed for YUZBEC
            ligand_a, ligand_b = sorted(ligands, key = lambda x: len(x), reverse = True)[:2]
        row = {'mol_id': mol_id, 'ligand_a': ligand_a, 'ligand_b': ligand_b}
        rows.append(row)
    
    if chemdraw:
        chemdraw_label = 'chemdraw'
    else:
        chemdraw_label = 'rdkit'

    df = pd.DataFrame(rows)
    # Sort it by the 'mol_id' column
    df = df.sort_values('mol_id')
    df.to_csv(f'ligands_{chemdraw_label}.csv', index = False)
