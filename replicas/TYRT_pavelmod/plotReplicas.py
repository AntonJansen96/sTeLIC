#!/usr/bin/env python3

import os
import matplotlib
import matplotlib.pyplot as plt
from science.parsing import loadxvg
import numpy as np

matplotlib.rcParams.update({'font.size': 24})

plt.figure(dpi=200)

folder = "replicas_3"

for num in range(1, 11):

    os.chdir(f"{folder}/s_{num}")

    # if (not os.path.exists("cphmd-coord-1.xvg")) or (not os.path.exists("MD.gro")):
    os.system("gmx cphmd -s MD.tpr -e MD.edr -numplot 1")

    data = loadxvg("cphmd-coord-1.xvg")[1]

    hist, bins = np.histogram(data, bins=35, range=(-0.10, 1.10), density=True)
    bins = [(bins[i] + bins[i + 1]) / 2 for i in range(0, len(bins) - 1)]
    plt.plot(bins, hist)

    os.chdir('../..')

plt.xlabel(r'$\lambda$-coordinate')
plt.ylabel('Probability density')
plt.tight_layout()
plt.savefig(f"{folder}.png")
