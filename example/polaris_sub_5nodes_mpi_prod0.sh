#!/bin/bash
#PBS -N nmnat2_openmm
#PBS -l select=5
#PBS -l walltime=01:00:00
#PBS -q debug-scaling
#PBS -l filesystems=home:grand:eagle
#PBS -A datascience
#PBS -o logs/
#PBS -e logs/
#PBS -m abe
#PBS -M avasan@anl.gov

module use /soft/modulefiles
module load conda
conda activate /lus/eagle/projects/datascience/avasan/envs/sst_tf216

PWD=/lus/eagle/projects/datascience/avasan/Simulations/NMNAT-2/Simulations_Dimer
cd $PWD

mpiexec -n 4 -ppn 4 --cpu-bind verbose,list:0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16:17,18,19,20,21,22,23,24,25,26,27,28,29,30,31:32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47:48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63 python run_openmm_mpi.py -T "prod0" -P ${PWD} > logs/openmm_mpi.prod0.log 2> logs/openmm_mpi.prod0.err

#for i in $(seq 1 5); do
#    NODE=$(sed -n "${i}p" $PBS_NODEFILE)
#    for j in $(seq 1 4); do
#        let rid=(${i}-1)*4+${j}
#        ssh $NODE "cd /lus/eagle/projects/datascience/avasan/Simulations/NMNAT-2/Simulations_Dimer/replica${rid}; python production.1.py" &
#    done
#done

# Wait for all background jobs to finish
wait

echo "All jobs are done!"
