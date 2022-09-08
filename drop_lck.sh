#!/bin/bash
#SBATCH -c 8
#
# number of nodes
#SBATCH -N 1
#
# node
#SBATCH -p defq
#
# set the $OMP_NUM_THREADS variable
ompthreads=$SLURM_JOB_CPUS_PER_NODE
export OMP_NUM_THREADS=$ompthreads

# working directory setup
mkdir /nobackup/b6019832/$1
cp config.json /nobackup/b6019832/$1
cp ./gp_lck /nobackup/b6019832/$1
mv ./slurm-${SLURM_JOB_ID}.out /nobackup/b6019832/$1
cd /nobackup/b6019832/$1

./gp_eqm
