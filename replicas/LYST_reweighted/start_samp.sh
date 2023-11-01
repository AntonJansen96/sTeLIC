#!/bin/bash
for (( i=2; i<3; i++ ))
do
    rm -r s_${i}
    mkdir s_${i}
    cd s_${i}
    cp ../md_cph.mdp .
    cp ../subm* .
    /projappl/project_2000996/gromacs-constantph/build/bin/gmx_mpi grompp -p ../topol.top -f md_cph.mdp -c ../8_npt.gro -n ../index.ndx -o run.tpr -maxwarn 1
    sbatch submit_puhti.sh
    cd ../
done

exit;

