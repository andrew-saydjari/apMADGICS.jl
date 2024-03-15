## This script takes a bunch of theoretical model spectra and creates samples of stellar lines for the starLines prior and then builds that prior using a Krylov decomposition.
## This is really a one node script, multi-node is not really necessary for this task.
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
    nsub_out = 50
    normPercent = 94 #nothing turns it off
    # normPercent = nothing #94 #nothing turns it off

    nsub_rnd1 = 60
    rnd1size = 800

    rnd_seed = 103
    th_LSF_R = 22_500

    # Prior Dictionary
    prior_dict = Dict{String,String}()
    prior_base = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/"
    # StarCont Samples
    prior_dict["korg_run_path"] = prior_base*"starLine_disk_KnackedKorg/"
    # prior_dict["out_dir"] = prior_base*"starLine_priors_norm94/"
    prior_dict["out_dir"] = prior_base*"starLine_priors/"

    prior_dict["LSF_mat_APO"] = prior_dir0*"2023_04_01/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir0*"2023_04_07/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_07 by AKS
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
    x_model = 15000:(1//100):17000;
end

@everywhere begin
    function grab_spectrum(fname;normPercent=nothing)
        if !isnothing(normPercent)
            dat = readdlm(fname)[:,2]
            return dat./percentile(dat,normPercent).-1
        else
            return readdlm(fname)[:,2].-1
        end
    end

    function samp_fxn_mult(matList,precompList,x)
        Vi = matList[1]
        return Vi*(Vi'*x)
    end

    function kryburyCompress_noDiag(M::LowRankMultMat,nfeat::Int;nsub=2,tol=1e-15)
        λ, V, info = eigsolve(M,
            nfeat,
            nsub,
            krylovdim=2*nsub,
            ishermitian=true,
            issymmetric=true,
            :LM,
            tol=tol
        );
        return hcat(V[1:nsub]...)*Diagonal(sqrt.(λ[1:nsub]))
    end

    function prelim_decomp(subsamples;nsub=60,normPercent=nothing)
        outsamp = zeros(length(x_model),length(subsamples))
        for (i, fval) in enumerate(subsamples)
            outsamp[:,i].=grab_spectrum(fval,normPercent=normPercent)
        end
        outsamp./=sqrt(length(subsamples));
        subMat = LowRankMultMat([outsamp],[],samp_fxn_mult);
        try
            return kryburyCompress_noDiag(subMat,length(x_model);nsub=nsub);
        catch
            return NaN
        end
    end

end

# Get all High Res Samples
rng = MersenneTwister(rnd_seed);
flist = shuffle(rng,sort(glob("*/*.dat",prior_dict["korg_run_path"])));
println("Number of Synthetic Spectra Samples: ", length(flist))

# Round 1 of Decomposition
iterpart = Iterators.partition(flist,rnd1size)
println("Number of $rnd1size Spectra Partitions in Round 1: ", length(iterpart))

@everywhere prelim_decomp_bind(subsamples) = prelim_decomp(subsamples;nsub=nsub_rnd1,normPercent=normPercent)
pout = @showprogress pmap(prelim_decomp_bind,iterpart);

global bind = 1
outsamp = zeros(length(x_model),length(iterpart)*nsub_rnd1);
for i=1:size(pout,1)
    outsamp[.. ,bind:(bind+nsub_rnd1-1)] .= pout[i]
    global bind += nsub_rnd1
end

serialize(prior_dict["korg_run_path"]*"korg_rnd1_kry_"*string(nsub_rnd1)*".jdat",outsamp)

# Final Krylov Decomposition
normfac = size(outsamp,2)/nsub_rnd1
outsamp./=sqrt(normfac);
subMat = LowRankMultMat([outsamp],[],samp_fxn_mult);

Vout = kryburyCompress_noDiag(subMat,length(x_model);nsub=nsub_out);

fname = prior_dict["out_dir"]*"APOGEE_stellar_kry_$(nsub_out)_fullres.h5"
dirName = splitdir(fname)[1]
if !ispath(dirName)
    mkpath(dirName)
end
h5write(fname,"Vmat",Vout)

# Read in FullRes StarLine Prior on Each Worker (prep for convolution)
@everywhere begin
    fname =  prior_dict["out_dir"]*"APOGEE_stellar_kry_$(nsub_out)_fullres.h5"
    Vout = h5read(fname,"Vmat")
end

# Convolve High Res Star Model with LSF
@everywhere begin
    offrng = -5//10:(1//10):(4//10)
    Vsubpix = zeros(length(wavetarg),nsub_out,length(offrng));
    function convolve_highRes_starModel(fibernum)
        fill!(Vsubpix,0)
        for (findx, fstep) in enumerate(offrng)
            Ksp = if fibernum>300
                deserialize(prior_dict["LSF_mat_LCO"]*"$findx"*"_"*lpad(fibernum-300,3,"0")*".jdat");
            else
                deserialize(prior_dict["LSF_mat_APO"]*"$findx"*"_"*lpad(fibernum,3,"0")*".jdat");
            end
            nvecLSF = dropdims(sum(Ksp,dims=2),dims=2);
            Vsubpix[:,:,findx] .= (Ksp*Vout)./nvecLSF;
        end
        lname =  prior_dict["out_dir"]*"APOGEE_stellar_kry_$(nsub_out)_subpix_f"*lpad(fibernum,3,"0")*".h5"
        h5write(lname,"Vmat",Vsubpix)
        return
    end
end

pout = @showprogress pmap(convolve_highRes_starModel,1:600);

## Make some reference LSFs
fsteprng = (5//10):(-1//10):(-4//10) #that should be left to right (for LSFs only)
sstep = 6.0e-6
@showprogress for (findx, fstep) in enumerate(fsteprng)
    nname = "ref_lsfs/sp_reflsf_norm_$findx.jdat"
    if !isfile(fname)
        wavetarg_new = 10 .^range(start=(4.179-125*sstep+fstep*sstep),step=sstep,length=8575+125);
        sp_lsf = instrument_lsf_sparse_matrix(x_model,wavetarg_new,th_LSF_R);
        
        dirNameN = splitdir(nname)[1]
        if !ispath(dirNameN)
            mkpath(dirNameN)
        end
        serialize(nname,sp_lsf)
    end
end

# Write out subpixel convolved starLine prior with theoretical LSF
Vsubpix = zeros(length(wavetarg),nsub_out,length(fsteprng));
for (findx, fstep) in enumerate(fsteprng)
    Ksp = deserialize("ref_lsfs/sp_reflsf_norm_$findx.jdat");
    nvecLSF = dropdims(sum(Ksp,dims=2),dims=2); # this was missing from this step until 2024_02_26
    Vsubpix[:,:,findx] .= (Ksp*Vout)./nvecLSF; 
end
fname =  prior_dict["out_dir"]*"APOGEE_stellar_kry_$(nsub_out)_subpix_th_$(th_LSF_R).h5"
h5write(fname,"Vmat",Vsubpix)
