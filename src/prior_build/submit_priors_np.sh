#!/bin/bash
#SBATCH --account=sdss-np
#SBATCH --partition=sdss-shared-np
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64

#SBATCH --mem=0 #requesting all of the memory on the node

#SBATCH --time=96:00:00
#SBATCH --job-name=apMADGICS_priors
#SBATCH --output=%x_%j.out
#SBATCH --err=%x_%j.err

#SBATCH --mail-type=ALL
#SBATCH --mail-user=7155301634@vtext.com
# ------------------------------------------------------------------------------

echo $SLURM_JOB_NODELIST

# export JULIA_NUM_GC_THREADS=1
# export JULIA_NUM_THREADS=1
# export JULIA_CPU_THREADS=1

# julia +1.10.0 sample_sky.jl # 3k core-h, 7.7h on 6 nodes, 10% CPU usage (2 corrupted skySpec_tellDiv_ files had to be manually rm-ed)
# julia +1.10.0 build_skyCont.jl # 672 core-h, 1.75h on 6 nodes, 100% cpu usage [OOM possible with Krylov]
julia +1.10.0 build_skyLines.jl # on 6 nodes, 100% cpu usage [never use Krylov]
# julia +1.10.0 sample_Tfun.jl
# julia +1.10.0 sample_starCont.jl
# julia +1.10.0 build_starCont.jl [OOM possible with Krylov]
# julia +1.10.0 sample_Korg.jl # 966.4 core-h, 2.5h on 6 nodes, 34.8 core-s/spec, 100% cpu usage
# julia +1.10.0 build_starLines.jl # 40 min on 1 node, 50% cpu usage
# julia +1.10.0 build_DIB.jl # 145 core-h, 2.3 h on 1 node, 100% cpu usage
# julia +1.10.0 build_starLines_dd.jl # try on 3 nodes? new parallelism strategy to try to mitigate time


# Clean up logs and Report Timing
formatted_time=$(printf '%dd %dh:%dm:%ds\n' $(($SECONDS/86400)) $(($SECONDS%86400/3600)) $(($SECONDS%3600/60)) $(($SECONDS%60)))
echo "Job completed in $formatted_time"

mkdir -p slurm_logs
mv ${SLURM_JOB_NAME}_${SLURM_JOBID}.out slurm_logs/${SLURM_JOB_NAME}_${SLURM_JOBID}.out
mv ${SLURM_JOB_NAME}_${SLURM_JOBID}.err slurm_logs/${SLURM_JOB_NAME}_${SLURM_JOBID}.err

# Note: OOM to be had from precompute H2O cross-sections for sample_Korg.jl