'''Combine separate ligands into full molecules'''

import sys
import os
import pandas as pd
from rdkit.Chem.rdmolfiles import MolToMolFile, MolFromSmiles
from modify_organometallics import ligate, fragment

def assemble(ligands):
    '''Assemble ligands into a molecule'''
    fragmented_smiles = sum([fragment(ligand) for ligand in ligands], [])
    fragmented = [MolFromSmiles(smiles) for smiles in fragmented_smiles]
    if any(fragmented_mol is None for fragmented_mol in fragmented):
        raise ValueError('Can\'t parse fragment smiles')
    return ligate(fragmented)

def clean(smiles):
    '''Clean up ligand smiles'''
    return smiles.replace('`', '').strip()

inpath = sys.argv[1]
outdir = sys.argv[2]

# Create the output directory if it isn't there already
if not os.path.exists(outdir):
    os.makedirs(outdir)

# Hand-annotated dataframe containing smiles strings of structures
database = pd.read_csv(inpath)

# Loop over rows, and construct molecules
ligand_colnames = ['SMILES 1', 'SMILES 2', 'SMILES 3']
full_structure_colname = 'SMILES full'
for index, row in database.iterrows():
    mol_id = row['entry']
    # Check if there is a full structure, and if there is use that
    if pd.notna(row[full_structure_colname]):
        molecule = MolFromSmiles(clean(row[full_structure_colname]))
        if molecule is None:
            print(f'Couldn\'t read molecule from molecule {mol_id}')
            continue
    # Otherwise, construct ligands from the ligand SMILES columns
    else:
        ligand_smiles = [row[colname] for colname in ligand_colnames if pd.notna(row[colname])]
        if len(ligand_smiles) == 0:
            continue
        if any('Ir' not in this_smiles for this_smiles in ligand_smiles):
            print(f'Missing iridium in molecule {mol_id}')
            continue
        ligands = [MolFromSmiles(clean(smiles)) for smiles in ligand_smiles]
        if any(ligand is None for ligand in ligands):
            print(f'Couldn\'t read ligand from molecule {mol_id}')
            continue
        try:
            molecule = assemble(ligands)
        except ValueError:
            print(f'Couldn\'t assemble ligands from {mol_id}')
            continue
    # Write molecule to the output directory
    outpath = os.path.join(outdir, f'{mol_id}.mol')
    MolToMolFile(molecule, outpath)
