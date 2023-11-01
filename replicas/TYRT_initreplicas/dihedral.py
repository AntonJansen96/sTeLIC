#!/usr/bin/env python3

import os
import matplotlib.pyplot as plt
import matplotlib

from science.utility import createIndexFile
from science.parsing import loadxvg

matplotlib.rcParams.update({'font.size': 18})
plt.figure(dpi=200)

# path = "s_1/MD.gro"
# createIndexFile(path, 'dihedral.ndx', ['resname TYRT and name CE2 CZ OH HH'])

for idx in range(1, 11):

    path = f"./s_{idx}/MD.xtc"
    os.system(f"gmx angle -f {path} -n CE1.ndx -of angdist.xvg -type dihedral -binwidth 10")

    data = loadxvg('angdist.xvg')
    plt.plot(data[0], data[1])

plt.xlabel("Degrees")
plt.ylabel("Probability density")
plt.tight_layout()
plt.savefig('CE1.png')
