#!/usr/bin/env python3

import numpy as np
import matplotlib
import matplotlib.pyplot as plt

from science.cphmd import BiasPotential
from science.cphmd import InverseBoltzmann
from science.parsing import loadxvg

coordsArray = []
for num in range(1, 11):
    coordsArray.append(loadxvg(f"./TYRT_initreplicas/s_{num}/cphmd-coord-1.xvg")[1])

obj = InverseBoltzmann('test', coordsArray, Nbins=120)
obj.plot('orig')
obj.addBiasPotential()
obj.plot('new')
# obj.addBiasPotential()
# obj.addpHPotential(pH=4.5, pKa=4.5)
# obj.plot('new')


# matplotlib.rcParams.update({'font.size': 24})

# bias = BiasPotential(dwpE=7.5)

# lamdas = np.arange(-0.08, 1.08, 0.01)
# energies = []
# for val in lamdas:
#     energies.append(bias.potential(val))

# plt.figure(dpi=200)
# plt.plot(lamdas, energies)
# plt.xlabel(r'$\lambda$-coordinate')
# plt.ylabel('Energy (kJ/mol)')
# plt.tight_layout()
# plt.savefig('barrier.png')
