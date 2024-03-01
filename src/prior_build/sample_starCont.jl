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
    using SortFilters, BasisFunctions, Random, DustExtinction, DelimitedFiles
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    runlist_range = 295 #1:600 #295, 245, 335, 101

    nsamp = 10_000

    Teff_rng = 4_000:1:10_000
    Av_rng = 0:1e-4:5
    Rv_rng = 2.6:1e-4:3.6

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Input List (not really a prior, but an input file we search for stars conditioned on)
    prior_dict["tellSamples2read"] = prior_dir*"2024_02_21/outlists/summary/dr17_dr17_tellSamples2read.txt"
    # Data for Cals (not really a prior, but an input the results depend on in detail)
    prior_dict["LSF_mat_APO"] = prior_dir0*"2023_04_01/mat_lsf_out/sp_combolsfmat_norm_6_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir0*"2023_04_07/mat_lsf_out/sp_combolsfmat_norm_6_" # last made 2023_04_07 by AKS
    prior_dict["fracTellSamples_APO"] = prior_dir0*"2023_04_03/outsamptell_apo.jdat" # last made 2023_04_03 by AKS
    prior_dict["fracTellSamples_LCO"] = prior_dir0*"2023_04_07/outsamptell_lco.jdat" # last made 2023_04_07 by AKS

    # Location of the Tfun samples
    tell_base = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/"
    prior_dict["tfun_samples"] = tell_base*"tell_prior_disk/tfun_samples_"
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
    x_model = 15000:0.01:17000

    tellsamplesTxt = readdlm(prior_dict["tellSamples2read"],',')
    tellsample_lst = map(x->x[x.!=""],eachrow(tellsamplesTxt))
end

@everywhere begin
    function genModSamp(intup,tellFracSamples,TfunSamples,Ksp,nvecLSF)
        Teff,Av,Rv,Tfunindx,Tfracindx = intup
        bbs = blackbody.(Ref(Teff), x_model*1e-8);
        rvec = redden_mult(x_model,Av,Rv);
        return tellFracSamples[:,Tfracindx].*TfunSamples[:,Tfunindx].*((Ksp*(rvec.*bbs))./nvecLSF);
    end
end

@everywhere begin
    function gen_starCont_samples(adjfibindx;loc_parallel=false,seed=203)

        savename = "tell_prior_disk/starCont_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            Ksp = if adjfibindx>300
                deserialize(prior_dict["LSF_mat_LCO"]*lpad(adjfibindx-300,3,"0")*".jdat");
            else
                deserialize(prior_dict["LSF_mat_APO"]*lpad(adjfibindx,3,"0")*".jdat");
            end

            nvecLSF = dropdims(sum(Ksp,dims=2),dims=2); # used only in starCont sample gen

            # Get tell obs to use from disk
            sampcatlst = []
            for fib_subindx in tellsample_lst[adjfibindx]
                push!(sampcatlst,deserialize(prior_dict["tfun_samples"]*lpad(fib_subindx,3,"0")*".jdat"))
            end
            TfunSamples = hcat(sampcatlst...);

            # Load the fraction telluric model samples (10k random from visits, stacked frames)
            tellFracSamples = if adjfibindx > 300
                deserialize(prior_dict["fracTellSamples_LCO"])
            else
                deserialize(prior_dict["fracTellSamples_APO"])
            end

            # Generate StarCont Samples
            rng = MersenneTwister(seed)
            # draw over parameter space (consider marginalizing over the dust exponent)
            Teff_lst = rand(rng,Teff_rng,nsamp)
            Av_lst = rand(rng,Av_rng,nsamp)
            Rv_lst = rand(rng,Rv_rng,nsamp)
            Tfunindx_lst = rand(rng,1:size(TfunSamples,2),nsamp);
            Tfracindx_lst = rand(rng,1:size(tellFracSamples,2),nsamp);
            itobj = Iterators.zip(Teff_lst,Av_lst,Rv_lst,Tfunindx_lst,Tfracindx_lst)

            genModSamp_bound(itobj) = genModSamp(itobj,tellFracSamples,TfunSamples,Ksp,nvecLSF)
            pout = if loc_parallel
                @showprogress pmap(genModSamp_bound,itobj); #not very fast because of passing
            else
                map(genModSamp_bound,itobj);
            end
            outsamp = zeros(length(wavetarg),size(pout,1));
            for i=1:size(pout,1)
                outsamp[:,i].=pout[i]
            end
            serialize(savename,outsamp)
        end
    end
end

# gen_starCont_samples(runlist_range,loc_parallel=true)
@showprogress pmap(gen_starCont_samples,1:600)