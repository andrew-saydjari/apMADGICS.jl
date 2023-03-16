#!/bin/bash
#SBATCH --account=sdss-kp
#SBATCH --partition=sdss-kp
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16

#SBATCH --mem-per-cpu=3750

#SBATCH --time=24:00:00
#SBATCH --job-name=apMADGICS
#SBATCH --output=%x_%j.out
#SBATCH --err=%x_%j.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=7155301634@vtext.com
# ------------------------------------------------------------------------------

module load julia

julia pipeline.jl
