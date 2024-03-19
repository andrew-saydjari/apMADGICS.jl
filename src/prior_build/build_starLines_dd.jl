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
println("Running Main on ", gethostname()); flush(stdout)

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
    using SortFilters, Random, KrylovKit, KryburyCompress, Glob, DelimitedFiles, DustExtinction
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext (# not clear if this causes issues in 1.10.2)
# slurm_cpu_lock()
# t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

println(BLAS.get_config()); flush(stdout)
using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    nsub_out = 50
    fracdatcut = 0.85
    datmskexpand = 20

    release_dir = "dr17"
    redux_ver = "dr17"
    dr_number = 17

    DRP_SNR_CUT = 80
    flux_conserve_cut = 1e-3
    snr_proxy_cut = 30
    sfd_cut = 0.05

    # See https://www.sdss4.org/dr17/irspec/apogee-bitmasks/
    bad_extratarg_bits = 2^1+2^2+2^3 #COMISSIONING (1), TELLURIC (2), APO1M (3)
    bad_starflag_bits = 2^0+2^1+2^3+2^4+2^12+2^13+2^16+2^17+2^18+2^19+2^20+2^21+2^22+2^23+2^25 
    #BAD_PIXELS (0), COMMISSIONING (1), VERY_BRIGHT_NEIGHBOR (3), LOW_SNR (4), PERSIST_JUMP_POS (12), PERSIST_JUMP_NEG (13), SUSPECT_RV_COMBINATION (16), SUSPECT_BROAD_LINES (17), BAD_RV_COMBINATION (18), RV_REJECT (19), RV_SUSPECT (20), MULTIPLE_SUSPECT (21), RV_FAIL (22), SUSPECT_ROTATION (23), MTPFLUX_LT_50 (25)
    bad_aspcap_bits = 2^0+2^1+2^3+2^4+2^5+2^8+2^10
    +2^11+2^12+2^13+2^16+2^17+2^18+2^19+2^20+2^21
    +2^22+2^23+2^24+2^26+2^27+2^28+2^29+2^30+2^31
    +2^32+2^33+2^34+2^35+2^36+2^40+2^41
    # TEFF_WARN (0), LOGG_WARN (1), M_H_WARN (3), ALPHA_M_WARN (4), C_M_WARN (5), CHI2_WARN (8), ROTATION_WARN (10), SN_WARN (11), SPEC_HOLE_WARN (12), ATMOS_HOLE_WARN (13), TEFF_BAD (16), LOGG_BAD (17), VMICRO_BAD (18), M_H_BAD (19), ALPHA_M_BAD (20), C_M_BAD (21), N_M_BAD (22), STAR_BAD (23), CHI2_BAD (24), ROTATION_BAD (26), SN_BAD (27), SPEC_HOLE_BAD (28), ATMOS_HOLE_BAD (29), VSINI_BAD (30), NO_ASPCAP_RESULT (31), MISSING_APSTAR (32), NO_GRID (33), BAD_FRAC_LOWSNR (34), BAD_FRAC_BADPIX (35), FERRE_FAIL (36), PROBLEM_TARGET (40), MULTIPLE_SUSPECT (41)  

    offrng = (-5//10):(1//10):(4//10)

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    prior_dict["past_run"] = prior_dir*"2024_03_16/outdir_wu_th/apMADGICS_out.h5" # update for the new run
    prior_dict["map2visit"] = prior_dir*"2024_03_05/outlists/summary/dr17_dr17_map2visit_1indx.h5"
    prior_dict["map2star"] = prior_dir*"2024_03_05/outlists/summary/dr17_dr17_map2star_1indx.h5"

    prior_dict["starLines_LSF"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/starLine_priors_norm94/APOGEE_stellar_kry_50_subpix_f"
    prior_dict["out_dir"] = prior_dir*"2024_03_16/apMADGICS.jl/src/prior_build/starLine_priors_norm94_dd/"
end

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg);
    x_model = 15000:(1//100):17000;

    offrng_flt = Float64.(offrng) # Lanczos interpolation wants the type to be Float, not Rational (probably my fault in Interpolations.jl)

    Vsubpix = zeros(length(wavetarg),nsub_out,length(offrng));
    pixmskMat = zeros(Bool,length(wavetarg),nsub_out,length(offrng));
end

@everywhere begin
    RV_pixoff_final = reader(prior_dict["past_run"],"RV_pixoff_final")

    keyval = "x_residuals_z_v0"
    savename_sub = chop(prior_dict["past_run"],tail=3)*"_"*keyval*".h5"
    f = h5open(savename_sub)
    keyvalres = "x_residuals_v0"
    savename_sub_res = chop(prior_dict["past_run"],tail=3)*"_"*keyvalres*".h5"
    g = h5open(savename_sub_res)
    keyvalref = "x_starContinuum_v0"
    savename_sub_ref = chop(prior_dict["past_run"],tail=3)*"_"*keyvalref*".h5"
    h = h5open(savename_sub_ref)
    keyvallineCof = "x_starLineCof_v0"
    savename_sub_starLineCof = chop(prior_dict["past_run"],tail=3)*"_"*keyvallineCof*".h5"
    m = h5open(savename_sub_starLineCof)

    function finalize_xnorm(intup,RV_pixoff_final,pf,V_starbasis)
        (indx, findx) = intup
        svald = RV_pixoff_final[findx]
        rval = indInt(svald)
        tval = indTenth(svald)
        return pf["x1normMat"][:,indx] .+ (V_starbasis*pf["x2normMat"][:,indx]).*pf["ynormMat"][:,indx]
    end

    function grab_star_spec(findx,RV_pixoff_final,f,g,h,m)
        svald = RV_pixoff_final[findx]
        
        sig = g[keyvalres][:,findx]./f[keyval][:,findx]
        subspec = replace((g[keyvalres][:,findx])./(sig.^2),NaN=>0)
        x1norm = shiftHelper(subspec,svald)
        
        x2norm_noshift = m[keyvallineCof][:,findx]
        
        subspec_ref = replace(h[keyvalref][:,findx]./(sig.^2),NaN=>0)
        ynorm = shiftHelper(subspec_ref,svald)
        return x1norm, x2norm_noshift, ynorm
    end
end

# for each fiber, build a ddmodel for the starLine component
@everywhere begin
    function solve_star_ddmodel_fiber(adjfibindx,clean_inds)

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

        finalize_xnorm_partial(indx2run) = finalize_xnorm(indx2run,RV_pixoff_final,pfhand,V_starbasis[:,:,6])

        if (!isfile(fname) & !isfile(fname_subpix))
            pout = map(finalize_xnorm_partial,enumerate(clean_inds));

            xnormMat = zeros(length(wavetarg),length(pout))
            for (subindx, subout) in enumerate(pout)
                xnormMat[:,subindx] .= subout
            end

            nodat = dropdims(sum(iszero.(ynormMat),dims=2),dims=2)./size(ynormMat,2);
            msk_starCor = expand_msk(nodat.<fracdatcut,rad=datmskexpand);

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

            for (sindx, shift) in enumerate(offrng_flt)
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
        close(pfhand)
    end
end

## Make clean inds and write out (add some flagging with our new metrics)

if !isfile(prior_dict["out_dir"]*"clean_inds.h5")
    # handle making and reading map2star/map2visit
    # (what to do if those things are not done?, would need to change cuts, later version question)
    map2star = h5read(prior_dict["map2star"],"map2star_1indx")
    map2visit = h5read(prior_dict["map2visit"],"map2visit_1indx");

    f = FITS(summary_file_by_dr(release_dir,redux_ver,dr_number,"allVisit"))
    GLON = read(f[2],"GLON")[map2visit]
    GLAT = read(f[2],"GLAT")[map2visit]
    VREL = read(f[2],"VREL")[map2visit];
    SNR = read(f[2],"SNR")[map2visit];
    STARFLAG = read(f[2],"STARFLAG")[map2visit]
    STARFLAGS = read(f[2],"STARFLAGS")[map2visit];
    close(f)

    f = FITS(summary_file_by_dr(release_dir,redux_ver,dr_number,"allStar"))
    EXTRATARG = read(f[2],"EXTRATARG")[map2star]
    ASPCAPFLAG = read(f[2],"ASPCAPFLAG")[map2star]
    FE_H_FLAG = read(f[2],"FE_H_FLAG")[map2star]

    TEFF = read(f[2],"TEFF")[map2star]
    LOGG = read(f[2],"LOGG")[map2star]
    X_H = read(f[2],"X_H")[map2star]
    M_H = read(f[2],"M_H")[map2star]
    VMICRO = read(f[2],"VMICRO")[map2star]
    VSINI = read(f[2],"VSINI")[map2star]
    ALPHA_M = read(f[2],"ALPHA_M")[map2star];
    close(f)

    ## This biases our DD model towards things that FERRE did well on, we should roll this back at some point
    ## To roll it back, we need to convince ourselves that we can handle the outliers well. v0.10
    ## I should visualize the Kiel diagram density changes with and without this cut.
    ## SNR similarly is a stellar type bias
    TARG_mask = .!(EXTRATARG .& bad_extratarg_bits.!=0);
    STARFLAG_masks = .!(STARFLAG .& bad_starflag_bits.!=0);
    ASPCAP_masks = .!(ASPCAPFLAG .& bad_aspcap_bits.!=0);
    Teff_masks = (.!isnan.(TEFF));
    FE_masks = (FE_H_FLAG .== 0);
    apg_msk = TARG_mask .& STARFLAG_masks .& ASPCAP_masks .& Teff_masks .& FE_masks
    println("Visits after TARG/STAR/ASPCAP/TEFF/FE masks: $(count(apg_msk)), $(100*count(apg_msk)/length(apg_msk))"); flush(stdout)

    ## Query SFD
    sfd_map = SFD98Map()
    sfd_reddening = sfd_map.(deg2rad.(GLON),deg2rad.(GLAT))
    sfd_msk = (sfd_reddening.<sfd_cut)

    ## MADGICS Cuts
    avg_flux_conservation = reader(prior_dict["past_run"],"avg_flux_conservation")
    adjfiberindx_vec = reader(prior_dict["past_run"],"adjfiberindx")
    msk_flux_conserve = (avg_flux_conservation .< flux_conserve_cut);
    RV_minchi2_final = reader(prior_dict["past_run"],"RV_minchi2_final")
    snr_proxy = sqrt.(-RV_minchi2_final);
    msk_MADGICS_snr = (snr_proxy .> snr_proxy_cut); 

    clean_msk = (adjfiberindx_vec.<=300) 
    clean_msk .&= (apg_msk .& (SNR.>DRP_SNR_CUT)) # Cuts on ASPCAP/Upstream Processing/Targetting
    clean_msk .&= sfd_msk # SFD Mask (low-reddening)
    clean_msk .&= (msk_flux_conserve .& msk_MADGICS_snr) # Remove StarConts that failed to converge well and low SNR model detections
    # clean_msk .&= (.!(-0.8 .< RV_pixoff_final .< 1.2)) .& ## add back moon avoidance only if we see deviation in the theory work up near zero
    println("Clean APO Visits for DD Model Training: $(count(clean_msk)), $(100*count(clean_msk)/count(adjfiberindx_vec.<=300))"); flush(stdout)
    clean_inds = findall(clean_msk);
    h5write(prior_dict["out_dir"]*"clean_inds.h5","clean_inds_apo",clean_inds)

    clean_msk = (adjfiberindx_vec.>300) 
    clean_msk .&= (apg_msk .& (SNR.>DRP_SNR_CUT)) # Cuts on ASPCAP/Upstream Processing/Targetting
    clean_msk .&= sfd_msk # SFD Mask (low-reddening)
    clean_msk .&= (msk_flux_conserve .& msk_MADGICS_snr) # Remove StarConts that failed to converge well and low SNR model detections
    # clean_msk .&= (.!(-0.8 .< RV_pixoff_final .< 1.2)) .& ## add back moon avoidance only if we see deviation in the theory work up near zero
    println("Clean LCO Visits for DD Model Training: $(count(clean_msk)), $(100*count(clean_msk)/count(adjfiberindx_vec.>300))"); flush(stdout)
    clean_inds = findall(clean_msk);
    h5write(prior_dict["out_dir"]*"clean_inds.h5","clean_inds_lco",clean_inds)
end

## make strip dd precursors and write out
@everywhere clean_inds_apo = h5read(prior_dict["out_dir"]*"clean_inds.h5","clean_inds_apo")
@everywhere clean_inds_lco = h5read(prior_dict["out_dir"]*"clean_inds.h5","clean_inds_lco")

@everywhere grab_star_spec_partial(indx2run) = grab_star_spec(indx2run,RV_pixoff_final,f,g,h,m)

# Run LCO
if !isfile(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5")
    pout = @showprogress pmap(grab_star_spec_partial,clean_inds_lco);

    x1normMat = zeros(length(wavetarg),length(pout))
    x2normMat = zeros(size(pout[1][2],1),length(pout))
    ynormMat = zeros(length(wavetarg),length(pout))
    for (subindx, subout) in enumerate(pout)
        x1normMat[:,subindx] .= subout[1]
        x2normMat[:,subindx] .= subout[2]
        ynormMat[:,subindx] .= subout[3]
    end

    h5write(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5","x1normMat",x1normMat)
    h5write(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5","x2normMat",x2normMat)
    h5write(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5","ynormMat",ynormMat)
    x1normMat = nothing
    x2normMat = nothing
    ynormMat = nothing
    GC.gc()
end

# Run APO
if !isfile(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5")
    pout = @showprogress pmap(grab_star_spec_partial,clean_inds_apo);

    x1normMat = zeros(length(wavetarg),length(pout))
    x2normMat = zeros(size(pout[1][2],1),length(pout))
    ynormMat = zeros(length(wavetarg),length(pout))
    for (subindx, subout) in enumerate(pout)
        x1normMat[:,subindx] .= subout[1]
        x2normMat[:,subindx] .= subout[2]
        ynormMat[:,subindx] .= subout[3]
    end

    h5write(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5","x1normMat",x1normMat)
    h5write(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5","x2normMat",x2normMat)
    h5write(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5","ynormMat",ynormMat)
    x1normMat = nothing
    x2normMat = nothing
    ynormMat = nothing
    GC.gc()
end

## Solve DD Model for StarLine Components (I split this part into 3 scripts 2 APO, 1 LCO because of how slow it was last time)

# Run LCO
@everywhere ynormMat = h5read(prior_dict["out_dir"]*"strip_dd_precursors_lco.h5","ynormMat")
GC.gc()
@everywhere solve_star_ddmodel_fiber_partial(adjfiberindx) = solve_star_ddmodel_fiber(adjfiberindx,clean_inds_lco)
@showprogress pmap(solve_star_ddmodel_fiber_partial,301:600) #obvi SVD speed up if we switch to MKL... do we really want two LA deps for this repo?

# Run APO (might need to write a workers per node handler script because we cannot load them all into memory at once)
@everywhere ynormMat = h5read(prior_dict["out_dir"]*"strip_dd_precursors_apo.h5","ynormMat")
GC.gc()
@everywhere solve_star_ddmodel_fiber_partial(adjfiberindx) = solve_star_ddmodel_fiber(adjfiberindx,clean_inds_apo)
@showprogress pmap(solve_star_ddmodel_fiber_partial,1:300) # tried switching the pmap outside, watch RAM

## check continuous connected components number in the msk_starCor
for tstind in 1:600
    zero_ranges = find_zero_ranges(V_subpix[:,1,1])
    println("Fiber $tstind has $(length(zero_ranges)) continuous connected zero components (should be 4).")
end