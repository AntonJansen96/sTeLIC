#!/bin/bash

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1
gmx mdrun -v -deffnm EM -npme 0

gmx grompp -f NVT.mdp -c EM.gro -p topol.top -n index.ndx -o NVT.tpr
gmx mdrun -v -deffnm NVT -npme 0

gmx grompp -f NPT.mdp -c NVT.gro -p topol.top -n index.ndx -o NPT.tpr
gmx mdrun -v -deffnm NPT -npme 0

# gmx grompp -f MD.mdp -c NPT.gro -p topol.top -n index.ndx -o MD.tpr
# gmx mdrun -v -deffnm MD -npme 0
