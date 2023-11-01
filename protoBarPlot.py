#!/usr/bin/env python3

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import pickle
import MDAnalysis
from science.utility import makeSuperDict
from science.parsing import loadxvg
from science.cphmd import protonation, deprotonation, getLambdaFileIndices

# PARAMETERS ###################################################################

u = MDAnalysis.Universe('/home/anton/stelic/production/closed_10/01/phneutral.pdb')
sims = ['closed_6.5', 'closed_10', 'open_6.5', 'open_10']
reps = [1, 2, 3]

residues = list(u.select_atoms('chainID A and resname ARGT LYST TYRT ASPT GLUT HSPT').residues.resids)
# print(residues)   # debug

# DATA PART ####################################################################

# We do this so we don't have to run the analysis part every time we plot.
if not os.path.exists('protonation_test.obj'):

    titratableHSPTids = list(u.select_atoms('chainID A and resname HSPT').residues.resids)
    # print(titratableHSPTids)  # debug

    protoData = makeSuperDict([sims, residues, []])
    # print(protoData)  # debug

    for target in residues:
        print(target)

        lambdaIndices = getLambdaFileIndices(u, target)
        # print(lambdaIndices)  # debug

        for sim in sims:
            for rep in reps:
                for idx in lambdaIndices:

                    x = loadxvg(f"/home/anton/stelic/production/{sim}/{rep:02d}/lambda_coords/cphmd-coord-{idx}.xvg", dt=100)[1]

                    # This is a bug fix. We should use the first out of three
                    # lambda_coordinates for HSPT as the first corresponds to single
                    # vs double protonated (see .mdp), however it is also reversed
                    # because of multistate! so use deprotonation instead:
                    if target in titratableHSPTids:
                        protoData[sim][target].append(deprotonation(x))
                    else:
                        protoData[sim][target].append(protonation(x))

    pickle.dump(protoData, open('protonation_test.obj', 'wb'))

# PLOTTING PART ################################################################

data = ['ASP', 'GLU', 'HIS', 'ARG', 'LYS', 'TYR'] + ['E7', 'D10', 'K16', 'D18', 'K26', 'E27', 'E28', 'R37', 'D39', 'R41', 'E48', 'H49', 'E53', 'K55', 'H56', 'R57', 'Y59', 'K66', 'E69', 'E70', 'K71', 'R74', 'Y80', 'H81', 'R86', 'D88', 'R92', 'E98', 'D99', 'Y104', 'E106', 'R107', 'D118', 'R120', 'D125', 'H132', 'D134', 'H140', 'R143', 'E146', 'D155', 'E159', 'E160', 'E161', 'E166', 'H170', 'H174', 'E176', 'K179', 'D181', 'R184', 'E188', 'H190', 'E192', 'R193', 'H194', 'Y197', 'Y198', 'R201', 'D221', 'Y222', 'K224', 'R225', 'D227', 'D246', 'R249', 'Y252', 'D257', 'R278', 'R279', 'E281', 'H283', 'K285', 'R290', 'K291', 'D293', 'Y295', 'Y300', 'Y304']

matplotlib.rcParams.update({'font.size': 24})

def stderr(array: list) -> float:
    return np.std(array) / np.sqrt(len(array))

# Load the means and standard errors from our pickle object.
protoMean = pickle.load(open('protonation_test.obj', 'rb'))
protoErr = makeSuperDict([sims, residues, 0])
for sim in sims:
    for residue in residues:
        protoErr[sim][residue] = stderr(protoMean[sim][residue])
        protoMean[sim][residue] = np.mean(protoMean[sim][residue])

# print(protoMean)  # debug
# print(protoErr)   # debug

ncols = 15
nrows = 6

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(1.8 * ncols, 4.5 * nrows), dpi=300)

row = 0
col = 0

for idx in range(0, len(data)):

    subplt = axs[row, col]

    width = 0.3
    x = np.arange(len([1]))

    if data[idx] in ['ASP', 'GLU', 'HIS', 'ARG', 'LYS', 'TYR']:
        mean2 = [{'ASP': 0.00, 'GLU': 0.00, 'HIS': 0.00, 'ARG': 1.00, 'LYS': 0.72, 'TYR': 0.41}[data[idx]]]
        mean1 = [{'ASP': 0.00, 'GLU': 0.01, 'HIS': 0.52, 'ARG': 1.00, 'LYS': 1.00, 'TYR': 1.00}[data[idx]]]

        subplt.bar(x - 0.5 * width, mean1, width, color='#43a2ca')
        subplt.text(x - 0.5 * width - 0.14, mean1[0] + 0.01, mean1[0], size=19)

        subplt.bar(x + 0.5 * width, mean2, width, color='#a8ddb5')
        subplt.text(x + 0.5 * width - 0.14, mean2[0] + 0.01, mean2[0], size=19)

    else:
        mean4 = [protoMean['closed_6.5'][float(data[idx][1:])]]
        serr4 = [protoErr['closed_6.5'][float(data[idx][1:])]]
        subplt.bar(     x - width * 1.5, mean4, width, color='#43a2ca')
        subplt.errorbar(x - width * 1.5, mean4, serr4, color='#43a2ca', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

        mean3 = [protoMean['closed_10'][float(data[idx][1:])]]
        serr3 = [protoErr['closed_10'][float(data[idx][1:])]]
        subplt.bar(     x - width / 2.0, mean3, width, color='#a8ddb5')
        subplt.errorbar(x - width / 2.0, mean3, serr3, color='#a8ddb5', fmt='none', capsize=7, linewidth=3, markeredgewidth=2)

        mean2 = [protoMean['open_6.5'][float(data[idx][1:])]]
        serr2 = [protoErr['open_6.5'][float(data[idx][1:])]]
        subplt.bar(     x + width / 2.0, mean2, width, color='#43a2ca', edgecolor='w', lw=0, hatch='//', zorder=13)
        subplt.errorbar(x + width / 2.0, mean2, serr2, color='#43a2ca', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=10)

        mean1 = [protoMean['open_10'][float(data[idx][1:])]]
        serr1 = [protoErr['open_10'][float(data[idx][1:])]]
        subplt.bar(     x + width * 1.5, mean1, width, color='#a8ddb5', edgecolor='w', lw=0, hatch='\\\\', zorder=12)
        subplt.errorbar(x + width * 1.5, mean1, serr1, color='#a8ddb5', fmt='none', capsize=7, linewidth=3, markeredgewidth=2, zorder=11)

    subplt.set_xticks([])           # Disable x-ticks and set x-label.
    subplt.set_xlabel(data[idx])
    subplt.set_ylim([0, 1.0])       # Set y-lim.

    # Remove boxes around plots.
    for loc in ['top', 'right', 'bottom', 'left']:
        subplt.spines[loc].set_visible(False)

    # If we're not in the first column, do not show the yticks.
    if col != 0:
        subplt.set_yticks([])
    else:
        subplt.set_ylabel('Protonation')
        subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])

    # Increment row and column indices.
    col += 1
    if col == ncols:
        row += 1
        col = 0

# Remove empty boxes if the last row wasn't filled.
if (nrows > 1) and (len(data) % ncols != 0):
    left = ncols - len(data) % ncols

    for ii in range(0, left):
        axs[nrows - 1, ncols - 1 - ii].axis('off')

fig.tight_layout(pad=0.8)
fig.savefig('stelic_barplot.png')
fig.clf()
