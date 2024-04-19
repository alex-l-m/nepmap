'''Search Cambridge Structure Database by doi read from stdin, and save result as a mol2 file'''

import sys
from os import mkdir
import ccdc.search
import ccdc.io

outdir = 'mol2'
try:
    mkdir(outdir)
except FileExistsError:
    pass

# Print header
print('entry,comments,doi,database,formal_charge,paper_url,database_url,deposition_number,molecule_formula,crystal_formula')
for doi in sys.stdin:
    doi = doi.strip()
    search = ccdc.search.TextNumericSearch()
    # Add doi to the search
    search.add_doi(doi)
    hits = search.search()
    for hit in hits:
        identifier = hit.identifier
        # The molecule is likely the largest component
        # Extract the molecule
        mol = hit.molecule.heaviest_component
        # Save in mol2 format
        outpath = f'{outdir}/{identifier}.mol2'
        ccdc.io.MoleculeWriter(outpath).write(mol)
        
        # Other properties I want for the output spreadsheet
        formal_charge = mol.formal_charge
        paper_url = f'https://doi.org/{doi}'
        database_url = f'https://www.ccdc.cam.ac.uk/structures/Search?Ccdcid={identifier}'
        deposition_number = hit.entry.ccdc_number
        molecule_formula = mol.formula
        crystal_formula = hit.crystal.formula

        # Output a row of data
        print(f'{identifier},,{doi},CSD,{formal_charge},{paper_url},{database_url},{deposition_number},{molecule_formula},{crystal_formula}')

