#!/usr/bin/env python3

import MDAnalysis

u = MDAnalysis.Universe("./charmm-gui-closed.pdb")

sel = list(u.select_atoms("resname ARG LYS TYR ASP GLU HSD").residues.resnames)

for res in ['ARG', 'LYS', 'TYR', 'ASP', 'GLU', 'HSD']:
    print(res, sel.count(res))

print(len(sel))
