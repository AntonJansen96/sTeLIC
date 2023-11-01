#!/usr/bin/env python3

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import MDAnalysis
from science.parsing import loadxvg

# PARAMETERS

structFile = "phneutral.pdb"  # This MUST have chain identifiers.
coordPath  = "../structures/EQ3/closed_6.5/HEAVY"
singlSite  = ['ASPT', 'GLUT']  + ['TYRT', 'ARGT', 'LYST']
multiSite  = ['HSPT']
numChains  = 5

################################################################################

matplotlib.rcParams.update({'font.size': 18})  # Matplotlib settings.

if not os.path.exists('boltzmann'):  # Create output dir if it doesn't exist.
    os.mkdir('boltzmann')

u = MDAnalysis.Universe(structFile)
sel = 'chainID A and resname ' + ' '.join(singlSite + multiSite)

# Lists of titratable resids and resnames in chain A.
resids = list(u.select_atoms(sel).residues.resids)
resnames = list(u.select_atoms(sel).residues.resnames)

# Total number of cphmd-coord files for chain A.
numCoordsA = len(resnames)
for name in multiSite:
    numCoordsA += 2 * resnames.count(name)

num = 1
for ii in range(0, len(resnames)):

    array = []
    for jj in range(0, numChains):
        array.append(num + jj * numCoordsA)

    if resnames[ii] in singlSite:
        num += 1

    elif resnames[ii] in multiSite:
        num += 3

    print(f"{resnames[ii]}-{resids[ii]}", array)

    data = []
    for coord in array:
        data += loadxvg(os.path.normpath(coordPath + f"/cphmd-coord-{coord}.xvg"))[1]

    plt.figure(dpi=200)
    hist, bins = np.histogram(data, bins=100, range=(-0.10, 1.10), density=True)
    bins = [(bins[i] + bins[i + 1]) / 2 for i in range(0, len(bins) - 1)]

    # This is the Boltzmann inversion part where we go from
    # probability density to energy U = - RT log(p):
    hist = [0.001 * -8.3145 * 300 * np.log(val) for val in hist]

    plt.plot(bins, hist)
    plt.hlines(y=min(hist), xmin=-0.1, xmax=1.1, linestyles='--', color='red')
    plt.grid()
    plt.title(f"Residue {resnames[ii]}-{resids[ii]} (all chains)")
    plt.xlim(-0.1, 1.1)
    plt.xlabel(r'$\lambda$-coordinate')
    plt.ylabel('Energy (kJ/mol)')
    plt.tight_layout()
    plt.savefig(f"boltzmann/{resids[ii]}.png")
    plt.close()
