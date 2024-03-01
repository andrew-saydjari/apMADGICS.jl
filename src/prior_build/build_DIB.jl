## This script builds the prior for DIBs using theoretical and analytic lineshapes, for now Gaussians.
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
    using SortFilters, Random, KrylovKit, KryburyCompress, Glob, DelimitedFiles
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    cwave_cent = 15273 # nominal center for DIB, pipeline.jl shifts prior to account for all other centers (assuming 15273 is the first DIB in the DIB list)

    amp = 3.0 # implemented in an EW-esque fashion
    damp = 0.2
    ddamp = 0.2

    # Prior Dictionary
    prior_dict = Dict{String,String}()
    prior_base = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/"

    prior_dict["out_dir"] = prior_base*"dib_priors/"

    prior_dict["LSF_mat_APO"] = prior_dir0*"2023_04_01/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir0*"2023_04_07/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_07 by AKS

    sigrng = 4//10:1//100:4//1
    offrng = -5//10:(1//10):(4//10)
end

@everywhere begin
    delLog = 6e-6; 
    wavetarg = 10 .^range(start=(4.179-125*delLog),step=delLog,length=8575+125)
    minw, maxw = extrema(wavetarg);
    x_model = 15000:(1//100):17000;
end

## Compute 1 and 3 component analytic models
Vall = zeros(length(wavetarg),3,length(sigrng),length(offrng))

for (sigindx,sigma) in enumerate(sigrng), (sindx,offset) in enumerate(offrng)
    x0 = cwave_cent*10^(offset*delLog)
    Vall[:,1,sigindx,sindx] .= gauss1d(amp/sqrt(sigma),x0,sigma,wavetarg)
    Vall[:,2,sigindx,sindx] .= gauss1d_deriv(damp*amp/sqrt(sigma),x0,sigma,wavetarg)
    Vall[:,3,sigindx,sindx] .= gauss1d_2deriv(ddamp*amp/sqrt(sigma),x0,sigma,wavetarg)
end

nkeep = 1
norm_array = dropdims(sum(Vall.^2,dims=1),dims=1);
covdet = zeros(length(sigrng),length(offrng))
for (sigindx,sigma) in enumerate(sigrng), (sindx,offset) in enumerate(offrng)
    covdet[sigindx,sindx] = logdet(2*pi.*Diagonal(norm_array[1:nkeep,sigindx,sindx]))
end

fname = prior_dict["out_dir"]*"precomp_dust_$(nkeep)_analyticDeriv_stiff.h5"
dirName = splitdir(fname)[1]
if !ispath(dirName)
    mkpath(dirName)
end
h5write(fname,"Vmat",Vall[:,1:nkeep,:,:])
h5write(fname,"covdet",covdet)

nkeep = 3
norm_array = dropdims(sum(Vall.^2,dims=1),dims=1);
covdet = zeros(length(sigrng),length(offrng))
for (sigindx,sigma) in enumerate(sigrng), (sindx,offset) in enumerate(offrng)
    covdet[sigindx,sindx] = logdet(2*pi.*Diagonal(norm_array[1:nkeep,sigindx,sindx]))
end

fname = prior_dict["out_dir"]*"precomp_dust_$(nkeep)_analyticDeriv_soft.h5"
h5write(fname,"Vmat",Vall[:,1:nkeep,:,:])
h5write(fname,"covdet",covdet) #but we wouldn't want to use this covdet

# Compute LSF convolved priors (for implicit deconvolution by projection in the latent space)
@everywhere begin
    Vall_lsf = zeros(length(wavetarg),3,length(sigrng),length(offrng))
    covdet_lsf = zeros(length(sigrng),length(offrng))
end

@everywhere begin
    function lsfDIB(adjfibindx)
        fill!(Vall_lsf,0)
        fill!(covdet_lsf,0)
        # analytic shift, so don't shift the LSF
        Ksp = if adjfibindx>300
            deserialize(prior_dict["LSF_mat_LCO"]*"6_"*lpad(adjfibindx-300,3,"0")*".jdat");
        else
            deserialize(prior_dict["LSF_mat_APO"]*"6_"*lpad(adjfibindx,3,"0")*".jdat");
        end
        nvecLSF = dropdims(sum(Ksp,dims=2),dims=2);

        for (sigindx,sigma) in enumerate(sigrng), (sindx,offset) in enumerate(offrng)
            x0 = cwave_cent*10^(offset*delLog)
            Vall_lsf[:,1,sigindx,sindx] .= (Ksp*gauss1d(amp/sqrt(sigma),x0,sigma,x_model))./nvecLSF
            Vall_lsf[:,2,sigindx,sindx] .= (Ksp*gauss1d_deriv(damp*amp/sqrt(sigma),x0,sigma,x_model))./nvecLSF
            Vall_lsf[:,3,sigindx,sindx] .= (Ksp*gauss1d_2deriv(ddamp*amp/sqrt(sigma),x0,sigma,x_model))./nvecLSF            
        end
        norm_array_lsf = dropdims(sum(Vall_lsf.^2,dims=1),dims=1);

        #
        nkeep=1
        for (sigindx,sigma) in enumerate(sigrng), (sindx,offset) in enumerate(offrng)
            covdet_lsf[sigindx,sindx] = logdet(2*pi.*Diagonal(norm_array_lsf[1:nkeep,sigindx,sindx]))
        end
        
        savename = prior_dict["out_dir"]*"precomp_dust_$(nkeep)_analyticDerivLSF_stiff_"*lpad(adjfibindx,3,"0")*".h5"
        h5write(savename,"Vmat",Vall_lsf[:,1:nkeep,:,:])
        h5write(savename,"covdet",covdet_lsf)
        
        #
        nkeep=3
        for (sigindx,sigma) in enumerate(sigrng), (sindx,offset) in enumerate(offrng)
            covdet_lsf[sigindx,sindx] = logdet(2*pi.*Diagonal(norm_array_lsf[1:nkeep,sigindx,sindx]))
        end
        
        savename = prior_dict["out_dir"]*"precomp_dust_$(nkeep)_analyticDerivLSF_soft_"*lpad(adjfibindx,3,"0")*".h5"
        h5write(savename,"Vmat",Vall_lsf[:,1:nkeep,:,:])
        h5write(savename,"covdet",covdet_lsf)
    end
end

pout = @showprogress pmap(lsfDIB,1:600);
