#!/bin/bash
#SBATCH -c 16
#
# number of nodes
#SBATCH -N 1
#
# set the $OMP_NUM_THREADS variable
ompthreads=$SLURM_JOB_CPUS_PER_NODE
export OMP_NUM_THREADS=$ompthreads

# working directory setup
mkdir /nobackup/b6019832/$1

#SBATCH –-workdir=/nobackup/b6019832/$1

cp ~/bb_qdrop_hom/config.json /nobackup/b6019832/$1
cp ~/bb_qdrop_hom/gp_eqm /nobackup/b6019832/$1

./gp_eqm

mv ~/bb_qdrop_hom/slurm-${SLURM_JOB_ID}.out /nobackup/b6019832/$1
