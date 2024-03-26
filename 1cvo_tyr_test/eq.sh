#!/bin/bash

source /usr/local/gromacs_constantph/bin/GMXRC

gmx grompp -f EM.mdp -c phneutral.pdb -p topol.top -n index.ndx -o EM.tpr -r phneutral.pdb -maxwarn 1
gmx mdrun -v -deffnm EM -c EM.pdb -npme 0

gmx grompp -f NVT.mdp -c EM.pdb -p topol.top -n index.ndx -o NVT.tpr -r EM.pdb
gmx mdrun -v -deffnm NVT -c NVT.pdb -npme 0

gmx grompp -f NPT.mdp -c NVT.pdb -p topol.top -n index.ndx -o NPT.tpr -r NVT.pdb
gmx mdrun -v -deffnm NPT -c NPT.pdb -npme 0
