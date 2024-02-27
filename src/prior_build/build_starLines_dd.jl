## This script takes a bunch of starLine components built using theoretical priors and creates a data-driven correction to the starLine prior.
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
    fracdatcut = 0.4

    offrng = -5//10:(1//10):(4//10)

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    prior_dict["past_run"] = prior_dir*"2023_09_26/outdir_wu_th/apMADGICS_out.h5" # update for the new run

    prior_dict["starLines_LSF"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/starLine_priors_norm94/APOGEE_stellar_kry_50_subpix_f"
    prior_dict["out_dir"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/starLine_priors_norm94_dd/"


@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
    x_model = 15000:0.01:17000;

    Vsubpix = zeros(length(wavetarg),nsub_out,length(offrng));
    pixmskMat = zeros(Bool,length(wavetarg),nsub_out,length(offrng));
end

@everywhere begin
    RV_pixoff_final = reader(prior_dict["past_run"],"RV_pixoff_final")

    function finalize_xnorm(intup,RV_pixoff_final,pf,V_starbasis)
        (indx, findx) = intup
        svald = RV_pixoff_final[findx]
        rval = indInt(svald)
        tval = indTenth(svald)
        return pf["x1normMat"][:,indx] .+ (V_starbasis*pf["x2normMat"][:,indx]).*pf["ynormMat"][:,indx]
    end
end

# for each fiber, build a ddmodel for the starLine component
@everywhere begin
    function solve_star_ddmodel_fiber(adjfibindx)

        f = h5open(prior_dict["starLines_LSF"]*lpad(adjfibindx,3,"0")*".h5")
        V_starbasis = read(f["Vmat"])
        close(f)

        fname = prior_dict["out_dir"]*"APOGEE_starCor_svd_"*string(nsub_out)*"_f"*lpad(adjfibindx,3,"0")*".h5"
        fname_subpix = prior_dict["out_dir"]*"APOGEE_starCor_svd_"*string(nsub_out)*"_subpix_f"*lpad(adjfibindx,3,"0")*".h5" 

        pfhand = if (adjfibindx .<=300)
            h5open(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5","r")
        else
            h5open(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5","r")
        end

        # could load in only if last fiber was not the same tele
        clean_inds = if (adjfibindx .<=300)
            h5read(prior_dict["out_dir"]*"clean_inds.h5","clean_inds_apo")
        else
            h5read(prior_dict["out_dir"]*"clean_inds.h5",,"clean_inds_lco")
        end;

        # could load in only if last fiber was not the same tele
        ynormMat = if (adjfibindx.<=300)
            h5read(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5","ynormMat")
        else
            h5read(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5","ynormMat")
        end;
        GC.gc()

        finalize_xnorm_partial(indx2run) = finalize_xnorm(indx2run,RV_pixoff_final,pfhand,V_starbasis[:,:,6])

        if (!isfile(fname) & !isfile(fname_subpix))
            pout = pmap(finalize_xnorm_partial,enumerate(clean_inds));

            xnormMat = zeros(length(wavetarg),length(pout))
            for (subindx, subout) in enumerate(pout)
                xnormMat[:,subindx] .= subout
            end

            nodat = dropdims(sum(iszero.(ynormMat),dims=2),dims=2)./size(ynormMat,2);
            msk_starCor = (nodat.<fracdatcut);

            xdat = xnormMat[msk_starCor,:]
            ydat = ynormMat[msk_starCor,:]
            Ccor = xdat*xdat'
            Ccor ./= (ydat*ydat');

            SF = svd(Ccor);
            EVEC = zeros(length(wavetarg),size(SF.U,2))
            EVEC[msk_starCor,:].=SF.U;

            ## Save results as prior (without subpix shifts)
            h5write(fname,"Vmat",EVEC[:,1:nsub_out]*Diagonal(sqrt.(SF.S[1:nsub_out])))
            h5write(fname,"λv",SF.S[1:nsub_out])
            h5write(fname,"msk_starCor",convert.(Int,msk_starCor))

            ## subpix shifts
            fill!(Vsubpix,0)
            fill!(pixmskMat,0)
            fltmsk = convert.(Float64,msk_starCor);

            for (sindx, shift) in enumerate(offrng)
                for eigind in 1:nsub_out
                    Vsubpix[:,eigind,sindx] .= (shiftHelper(EVEC[:,eigind],-shift)).*sqrt.(SF.S[eigind]);
                    shiftmsk = shiftHelper(fltmsk,-shift);
                    pixmskMat[:,eigind,sindx] .= (abs.(shiftmsk.-1).<1e-2)
                end
            end

            final_msk = dropdims(minimum(dropdims(minimum(pixmskMat,dims=3),dims=3),dims=2),dims=2)

            for (sindx, shift) in enumerate(offrng)
                for eigind in 1:nsub_out
                    Vsubpix[.!final_msk,eigind,sindx] .= 0;
                end
            end

            h5write(fname_subpix,"Vmat",Vsubpix)
            h5write(fname_subpix,"msk_starCor",convert.(Int,final_msk))
        end
    end
end

## make strip dd precursors and write out
## make clean inds and write out (add some flagging with our new metrics)

map(solve_star_ddmodel_fiber,1:600) #pmap inside, SVD speed up if we switch to MKL... do we really want two BLAS deps for this repo?

## check continuous connected components number in the msk_starCor (really only need to check 1 APO and 1 LCO)