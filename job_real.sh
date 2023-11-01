#!/bin/bash
#SBATCH --time=1-23:59:59
#SBATCH --nodes=1
#SBATCH -p lindahl3,lindahl4
#SBATCH --job-name=NAME
#SBATCH --mail-user=iday@kth.se
#SBATCH --mail-type=ALL
#SBATCH -G 1

module load gromacs/2021-cphmd-beta1

# Create the .tpr file if it doesn't exist.
if [ ! -e MD.tpr ]
then
    gmx grompp -f MD.mdp -c CA.gro -p topol.top -n index.ndx -o MD.tpr -maxwarn 1
fi

# Perform mdrun (for 0.99 * 48h, -maxh option).
gmx mdrun -deffnm MD -cpi MD.cpt -npme 0 -maxh 48 -bonded gpu -nt $SLURM_JOB_CPUS_PER_NODE

# Resubmit this jobscript if MD.pdb (last simulation frame) doesn't exist.
if [ ! -e MD.gro ]
then
    sbatch job_real.sh
fi
