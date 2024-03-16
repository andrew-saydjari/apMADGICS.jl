## This is script grabs a bunch of HOT STD star spectra and decomposes them into continuum and skyline components, to serve fit samples of the transfer function required for building the starContinuum prior
# Author - Andrew Saydjari, CfA

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("../../"); Pkg.instantiate(); Pkg.precompile()
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
using BLISBLAS
using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager(),exeflags=["--project=../../"])
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker allocation took $dt"); t_then = t_now; flush(stdout)

@everywhere begin
    using BLISBLAS
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using ThreadPinning

    prior_dir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/"
    prior_dir0 = "/uufs/chpc.utah.edu/common/home/u6039752/scratch/working/"
    src_dir = "../../"
    include(src_dir*"src/utils.jl")
    include(src_dir*"src/gridSearch.jl")
    include(src_dir*"src/componentAndPosteriors.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/ingest.jl")
    include(src_dir*"src/lowRankPrescription.jl")
    include(src_dir*"src/marginalizeEW.jl")
    include(src_dir*"src/spectraInterpolation.jl")
    include(src_dir*"src/chi2Wrappers.jl")
    include(src_dir*"src/prior_build/prior_utils.jl")
    
    using StatsBase, ProgressMeter
    using SortFilters, BasisFunctions, Random, DustExtinction
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext (# not clear if this causes issues in 1.10.2)
# slurm_cpu_lock()
# t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

println(BLAS.get_config()); flush(stdout)
using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    runlist_range = 295 #1:600 #295, 245, 335, 101

    nsub = 60

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # StarCont Samples
    star_base = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/"
    prior_dict["starcont"] = star_base*"tell_prior_disk/starCont_"

    sky_base = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/"
    prior_dict["chebmsk_exp"] = sky_base*"sky_prior_disk/chebmsk_exp_"
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
end

@everywhere begin    
    function build_starCont(adjfibindx)
        fname = "star_priors/APOGEE_starcont_svd_"*string(nsub)*"_f"*lpad(adjfibindx,3,"0")*".h5"
        if !isfile(fname)
            savename = prior_dict["starcont"]*lpad(adjfibindx,3,"0")*".jdat"
            starcont = deserialize(savename)

            savename = prior_dict["chebmsk_exp"]*lpad(adjfibindx,3,"0")*".jdat"
            chebmsk_exp = deserialize(savename);

            specsum = dropdims(sum(starcont,dims=1),dims=1)
            Vred = starcont[chebmsk_exp,specsum.>0];
            mnorm = mean(filter(.!iszero,Vred))
            Vred./=mnorm
            # weights = ones(size(Vred,2));
            # Vred .*= reshape(weights,1,:);
            nsamp = size(Vred,2)
            # norm_weights = weights'*weights
            Cexp = Vred*Vred'
            # Csky./=norm_weights
            Cexp./=nsamp

            SF = svd(Cexp);
            EVEC = zeros(length(wavetarg),size(SF.U,2))
            EVEC[chebmsk_exp,:].=SF.U;

            dirName = splitdir(fname)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            h5write(fname,"Vmat",EVEC[:,1:nsub]*Diagonal(sqrt.(SF.S[1:nsub])))
            h5write(fname,"λv",SF.S[1:nsub])
            h5write(fname,"chebmsk_exp",chebmsk_exp)
        end
    end
end

# build_starCont(runlist_range)
@showprogress pmap(build_starCont,1:600) # last run was 4.5h on 6 np nodes