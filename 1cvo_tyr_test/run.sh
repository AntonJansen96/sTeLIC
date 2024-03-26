#!/bin/bash

source /usr/local/gromacs_constantph/bin/GMXRC

# gmx grompp -f MD.mdp -c NPT.pdb -p topol.top -n index.ndx -o MD.tpr
gmx mdrun -v -deffnm MD -npme 0 -bonded gpu
