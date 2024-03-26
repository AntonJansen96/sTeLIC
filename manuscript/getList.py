#!/usr/bin/env python3

import MDAnalysis

u = MDAnalysis.Universe("/home/anton/stelic/production/closed_10/01/CA_chains.pdb")
sel1 = list(u.select_atoms("resname ASPT GLUT HSPT TYRT ARGT LYST and chainID A").residues.resids)
# sel2 = list(u.select_atoms("resname ASPT GLUT HSPT TYRT ARGT LYST and chainID A").residues.resnames)

for idx in range(len(sel1)):
    print(sel1[idx])
