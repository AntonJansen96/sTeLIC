#!/usr/bin/env python3

import os
import numpy as np
import MDAnalysis
import matplotlib
import matplotlib.pyplot as plt

from science.cphmd import getLambdaFileIndices, movingDeprotonation
from science.parsing import loadxvg, pickleDump, pickleLoad
from science.utility import makeSuperDict

matplotlib.rcParams.update({'font.size': 20})

def getFullResName(num: int) -> str:
    data = ['E7', 'D10', 'K16', 'D18', 'K26', 'E27', 'E28', 'R37', 'D39', 'R41', 'E48', 'H49', 'E53', 'K55', 'H56', 'R57', 'Y59', 'K66', 'E69', 'E70', 'K71', 'R74', 'Y80', 'H81', 'R86', 'D88', 'R92', 'E98', 'D99', 'Y104', 'E106', 'R107', 'D118', 'R120', 'D125', 'H132', 'D134', 'H140', 'R143', 'E146', 'D155', 'E159', 'E160', 'E161', 'E166', 'H170', 'H174', 'E176', 'K179', 'D181', 'R184', 'E188', 'H190', 'E192', 'R193', 'H194', 'Y197', 'Y198', 'R201', 'D221', 'Y222', 'K224', 'R225', 'D227', 'D246', 'R249', 'Y252', 'D257', 'R278', 'R279', 'E281', 'H283', 'K285', 'R290', 'K291', 'D293', 'Y295', 'Y300', 'Y304']
    for idx in range(0, len(data)):
        if num == int(data[idx][1:]):
            return data[idx]

def stderr(array):
    return np.std(array) / np.sqrt(len(array))

# resids = [7, 10, 16, 18, 26, 27, 28, 37, 39, 41, 48, 49, 53, 55, 56, 57, 59, 66, 69, 70, 71, 74, 80, 81, 86, 88, 92, 98, 99, 104, 106, 107, 118, 120, 125, 132, 134, 140, 143, 146, 155, 159, 160, 161, 166, 170, 174, 176, 179, 181, 184, 188, 190, 192, 193, 194, 197, 198, 201, 221, 222, 224, 225, 227, 246, 249, 252, 257, 278, 279, 281, 283, 285, 290, 291, 293, 295, 300, 304, 304]
# resids = [7, 10, 16, 18, 26, 27, 28, 37, 39, 41, 48, 49, 53, 55, 56, 57, 59, 66, 69, 70, 71, 74, 80, 81, 86, 88, 92, 98, 99, 104, 106, 107, 118, 120, 125, 132, 134, 140, 143, 146]
resids = [155, 159, 160, 161, 166, 170, 174, 176, 179, 181, 184, 188, 190, 192, 193, 194, 197, 198, 201, 221, 222, 224, 225, 227, 246, 249, 252, 257, 278, 279, 281, 283, 285, 290, 291, 293, 295, 300, 304, 304]
sims   = ['closed_6.5', 'closed_10', 'open_6.5', 'open_10']
fancy  = ["Closed, pH 6.5", "Closed, pH 10", "Open, pH 6.5", "Open, pH 10"]
reps   = [1, 2, 3]
chains = ['A', 'B', 'C', 'D', 'E']

# GATHER ALL DATA INTO ONER SUPERDICT

if os.path.isfile("protoTimeCompact.obj"):
    
    superData = pickleLoad("protoTimeCompact.obj")
    print("Loaded pickled protoTimeCompact object...")

else:

    u = MDAnalysis.Universe("/home/anton/stelic/production/closed_10/01/CA_chains.pdb")

    # Holds moving deprotonations for all resids, sims, reps, chains.
    superData = makeSuperDict([resids, sims, reps, chains, []])  

    for resid in resids:
        for sim in sims:
            print(resid, sim)
            for rep in reps:
                for chain in chains:
                    num = getLambdaFileIndices(u, resid)[chains.index(chain)]
                    data = loadxvg(f"/home/anton/stelic/production/{sim}/{rep:02d}/lambda_coords/cphmd-coord-{num}.xvg", dt=100)

                    if resid in [49, 56, 81, 132, 140, 170, 174, 190, 194, 283]:    # HIS
                        array = [1 - val for val in movingDeprotonation(data[1])]
                    else:                                                           # Rest
                        array = movingDeprotonation(data[1])

                    superData[resid][sim][rep][chain] = array

    pickleDump(superData, "protoTimeCompact.obj")

# GET THE AVERAGE OVER THE RUNNING AVERAGES OF THE FIVE CHAINS

superMean = makeSuperDict([resids, sims, reps, []])

for resid in resids:
    for sim in sims:
        for rep in reps:

            length = len(superData[resid][sim][rep]['A'])
            comb = [0] * length

            for idx in range(length):
                a = 1 - superData[resid][sim][rep]['A'][idx]  # 1 - x because we go 
                b = 1 - superData[resid][sim][rep]['B'][idx]  # from deprotonation
                c = 1 - superData[resid][sim][rep]['C'][idx]  # to protonation!
                d = 1 - superData[resid][sim][rep]['D'][idx]
                e = 1 - superData[resid][sim][rep]['E'][idx]

                comb[idx] = np.mean([a, b, c, d, e])

            superMean[resid][sim][rep] = comb

# GET THE AVERAGE OVER THE AVERAGE OF THE FOUR REPLICATES

ultraMean = makeSuperDict([resids, sims, []])
ultraSdev = makeSuperDict([resids, sims, []])

for resid in resids:
    for sim in sims:
        
        l1 = len(superMean[resid][sim][1])
        l2 = len(superMean[resid][sim][2])
        l3 = len(superMean[resid][sim][3])

        length = min([l1, l2, l3])
        
        comb = length * [0]
        sdev = length * [0]

        for idx in range(length):
            
            data = []
            for rep in reps:
                for chain in chains:
                    data.append(1 - superData[resid][sim][rep][chain][idx])

            comb[idx] = np.mean(data)
            sdev[idx] = stderr(data)

        ultraMean[resid][sim] = comb
        ultraSdev[resid][sim] = sdev

nrows = 10
ncols = 4

fig, axs = plt.subplots(nrows=nrows, ncols=ncols, figsize=(17, 23), dpi=200)

row = 0
col = 0

for resid in resids:

    subplt = axs[row, col]

    # Plot the mean of twenty chains
    for sim in sims:
        x = ultraMean[resid][sim]
        t = range(len(x))

        lower = [0] * len(t)
        upper = [0] * len(t)

        for idx in t:
            lower[idx] = x[idx] - ultraSdev[resid][sim][idx]
            upper[idx] = x[idx] + ultraSdev[resid][sim][idx]

        if sim == 'closed_6.5':
            subplt.plot(t, x, linewidth=2, color='#43a2ca')
            subplt.fill_between(t, lower, upper, color='#43a2ca', alpha=0.5)

        if sim == 'closed_10':
            subplt.plot(t, x, linewidth=2, color='#a8ddb5')
            subplt.fill_between(t, lower, upper, color='#a8ddb5', alpha=0.5)

        if sim == 'open_6.5':
            subplt.plot(t, x, linewidth=2, color='#43a2ca', linestyle='--')
            subplt.fill_between(t, lower, upper, color='#43a2ca', alpha=0.5, hatch='//', edgecolor='w', lw=1)
    
        if sim == 'open_10':
            subplt.plot(t, x, linewidth=2, color='#a8ddb5', linestyle='--')
            subplt.fill_between(t, lower, upper, color='#a8ddb5', alpha=0.5, hatch='\\\\', edgecolor='w', lw=1)

    # Set x-lim and y-lim.
    subplt.set_ylim([-0.05, 1.05])
    subplt.set_xlim([0, 2000])

    # If we're not in the last row, do not show the xticks.
    subplt.set_xticks([400, 800, 1200, 1600])
    if row != nrows - 1:
        subplt.set_xticklabels([])
        subplt.xaxis.set_ticks_position('none')
    else:
        subplt.set_xlabel("Time (ns)")

    # If we're not in the first column, do not show the yticks.
    subplt.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    if col != 0:
        subplt.set_yticklabels([])
        subplt.yaxis.set_ticks_position('none')
    else:
        subplt.set_ylabel("Proto frac")

    # Add (horizontal) grid to all subplots.
    subplt.grid(True, linestyle='--', axis='y')

    # Add title to each subplot.
    subplt.set_title(getFullResName(resid), size=20)

    # Increment row and column indices.
    col += 1
    if col == ncols:
        row += 1
        col = 0

fig.tight_layout(pad=0.2)
fig.savefig(f"allinone2.png")
os.system(f"convert allinone2.png -trim allinone2.png")
fig.clear()
