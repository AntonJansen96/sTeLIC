#/usr/bin/env python3

import MDAnalysis
import MDAnalysis.analysis.rms
import matplotlib
import matplotlib.pyplot as plt
import os

from science.parsing import loadxvg
from science.utility import gromacs
from science.utility import createIndexFile

matplotlib.rcParams.update({'font.size': 22})

basePath = "/home/anton/stelic/production/closed_10/01"

# GET CA.pdb, including all chains in working DIR
os.system(f"addChainIdentifiers.py {basePath}/CA.gro ./CA_chains.pdb")

#! Do the RMSD analysis using gmx rms

chains = ['A', 'B', 'C', 'D', 'E']

for chain in chains:
    createIndexFile("CA_chains.pdb", "rmsd.ndx", groups=[f"chainID {chain} and resid 1 to 196 and name CA"])
    gromacs(f"rms -f {basePath}/MD_conv.xtc -s {basePath}/MD.tpr -n rmsd.ndx -o rmsd.xvg")

    data = loadxvg("rmsd.xvg")
    x = [val / 1000.0 for val in data[0]]
    t = [val * 10 for val in data[1]]

    plt.plot(x, t, label=chain)

plt.legend()
plt.xlabel("Time (ns)")
plt.ylabel("RMSD (A)")
plt.axis([0, 1000, 0, 4])
plt.tight_layout()
plt.savefig('rmsdtest1.png')

plt.close()
plt.clf()

# Do the RMSD analysis using MDAnalysis

# u = MDAnalysis.Universe("CA_chains.pdb", f"{basePath}/MD_conv.xtc")

# for chain in chains:
#     R = MDAnalysis.analysis.rms.RMSD(u, select=f"chainID {chain} and resid 1 to 196")
#     R.run()

#     t = [val / 1000.0 for val in R.rmsd.T[1]]
#     x = R.rmsd.T[2]

#     plt.plot(t, x, label=chain)

# plt.legend()
# plt.xlabel("Time (ns)")
# plt.ylabel("RMSD (A)")
# plt.axis([0, 1000, 0, 4])
# plt.tight_layout()
# plt.savefig('rmsdtest2.png')
