#!/usr/bin/env python3

from science.parsing import Structure

sourceFile = "input.pdb"
outputFile = "output.pdb"

letters = ['E', 'D', 'C', 'B', 'A']

pdb = Structure(sourceFile)

ii = 0
try:
    for idx in range(0, len(pdb.d_residues)):
        if pdb.d_residues[idx].d_resname not in ['SOL', 'BUF', 'POPC', 'NA', 'CL', 'SOD', 'TIP3', 'CLA']:
            pdb.d_residues[idx].d_chain = letters[ii]

        if pdb.d_residues[idx].d_resid + 1 != pdb.d_residues[idx + 1].d_resid and ii != len(letters):
            ii += 1

except IndexError:
    pass

pdb.write(outputFile)
