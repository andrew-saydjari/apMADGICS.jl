#!/bin/bash
#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64

#SBATCH --mem-per-cpu=7500

#SBATCH --time=96:00:00
#SBATCH --job-name=apMADGICS
#SBATCH --output=%x_%j.out
#SBATCH --err=%x_%j.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=7155301634@vtext.com
# ------------------------------------------------------------------------------

module load julia

julia pipeline.jl

mkdir -p slurm_logs
mv ${SLURM_JOB_NAME}_${SLURM_JOBID}.out slurm_logs/${SLURM_JOB_NAME}_${SLURM_JOBID}.out
mv ${SLURM_JOB_NAME}_${SLURM_JOBID}.err slurm_logs/${SLURM_JOB_NAME}_${SLURM_JOBID}.err