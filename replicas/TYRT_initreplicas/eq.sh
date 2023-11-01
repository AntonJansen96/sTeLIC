#!/bin/bash

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -maxwarn 1
gmx mdrun -v -deffnm EM -npme 0

gmx grompp -f NVT.mdp -c EM.gro -p topol.top -n index.ndx -o NVT.tpr
gmx mdrun -v -deffnm NVT -npme 0

gmx grompp -f NPT.mdp -c NVT.gro -p topol.top -n index.ndx -o NPT.tpr
gmx mdrun -v -deffnm NPT -npme 0

for i in {1..10}; do mkdir s_${i}; cp -r charmm36-mar2019-cphmd.ff index.ndx jobscript.sh MD.mdp NPT.gro topol.top s_${i}; done
