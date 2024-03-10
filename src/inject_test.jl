## This generates a series of injection tests using synthetic stars an DIBs in real sky fiber/continuum observations in order to test and calibrate the apMADGICS pipeline.
# Author - Andrew Saydjari, CfA

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("../"); Pkg.instantiate(); Pkg.precompile()
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
using BLISBLAS
using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager(),exeflags=["--project=../"])
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
    src_dir = "../"
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
    using SortFilters, Random, Distributions, Glob, DelimitedFiles, PoissonRandom
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    nsamp = 5_080

    adjfibindx = 295 # sky samples injections are made into, running simulated observed on this fiber
    fiberindx = if adjfibindx>300
        adjfibindx-300
    else
        adjfibindx
    end

    RV_range_pix = (-68,68) # pixscale is ~4.14 km/s per pixel

    skycont_only = false
    no_sky = false
    dibs_on = true

    dib_center_lambda_lst = [15273] #,15672]
    dib_ew_range = (-1.5,0)
    dib_sig_range = (0.7,3.7)
    dib_vel_range = (-450, 450) # km/s

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    prior_dict["out_dir"] = prior_dir*"2024_03_08/inject_15273only_295_real/"
    prior_dict["inject_cache_dir"] = prior_dir*"2024_03_08/inject_local_cache_15273only_real/"
    prior_dict["local_cache"] = prior_dir*"2024_03_08/local_cache_inject_real/"

    prior_dict["past_run"] = prior_dir*"2024_03_08/outdir_wu_295_LocMean/apMADGICS_out.h5" # use component decomp to only inject into star component
    prior_dict["past_runlst"] = prior_dir*"2024_01_19/outlists/dr17_dr17_star_input_lst_msked_"
    prior_dict["map2visit"] = prior_dir*"2024_03_05/outlists/summary/dr17_dr17_map2visit_1indx.h5"
    prior_dict["map2star"] = prior_dir*"2024_03_05/outlists/summary/dr17_dr17_map2star_1indx.h5"
    dr_number = 17

    prior_dict["LSF_mat_APO"] = prior_dir0*"2023_04_01/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir0*"2023_04_07/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_07 by AKS

    rnd_seed = 695

    DRP_SNR_CUT = 35
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
end

println("NumSamples Acceptable for Filenaming: ",length(string(nsamp)) <= 7)

@everywhere begin
    delLog = 6e-6;
    wavetarg = 10 .^range(start = (4.179-125*delLog),step=delLog,length=8575+125)
    minw, maxw = extrema(wavetarg)
    x_model = 15000:(1//100):17000;
    c = 299792.458; # in km/s

    Xd_stack = zeros(3*2048)
    Xd_std_stack = zeros(3*2048)
    waveobs_stack = zeros(3*2048)
    waveobs_stack_old = zeros(3*2048)
    pixmsk_stack = zeros(Int,3*2048)
    telluric_stack = zeros(3*2048)
    fullBit = zeros(Int,3*2048);
    outvec = zeros(length(wavetarg))
    outvar = zeros(length(wavetarg))
    cntvec = zeros(Int,length(wavetarg));
end

@everywhere begin
    err_correct_Dict = deserialize(prior_dict["chip_fluxdep_err_correction"])

    ntuplst = deserialize(prior_dict["past_runlst"]*lpad(adjfibindx,3,"0")*".jdat")
    # these are semi gross globals; they are used in the take_draw function
    x_starLines = reader(prior_dict["past_run"],"x_starLines_v0");
    x_starContinuum = reader(prior_dict["past_run"],"x_starContinuum_v0");

    map2star_1indx = h5read(prior_dict["map2star"],"map2star_1indx")
    map2star_adjfibindx = h5read(prior_dict["map2star"],"adjfibindx")
    map2star = map2star_1indx[map2star_adjfibindx.==adjfibindx];
    
    map2visit = h5read(prior_dict["map2visit"],"map2visit_1indx");
    
    map2visit_1indx = h5read(prior_dict["map2visit"],"map2visit_1indx")
    map2visit_adjfibindx = h5read(prior_dict["map2visit"],"adjfibindx")
    map2visit = map2visit_1indx[map2visit_adjfibindx.==adjfibindx];
    
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
    println("Visits after TARG/STAR/ASPCAP/TEFF/FE masks: $(count(apg_msk)), $(100*count(apg_msk)/length(apg_msk))")
    
    ## Query SFD
    sfd_map = SFD98Map()
    sfd_reddening = sfd_map.(deg2rad.(GLON),deg2rad.(GLAT))
    sfd_msk = (sfd_reddening.<sfd_cut);

    clean_msk = (apg_msk .& (SNR.>DRP_SNR_CUT)) # Cuts on ASPCAP/Upstream Processing/Targetting
    clean_msk .&= sfd_msk # SFD Mask (low-reddening)
    println("Clean Visits for Injection Tests: $(count(clean_msk)), $(100*count(clean_msk)/length(clean_msk))")

    obs_indices2use = findall(clean_msk)
end

@everywhere begin
    function take_draw(ovtup; caching=true, dib_center_lambda_lst=[15273], inject_cache_dir="./inject_local_cache",cache_dir="./local_cache")
        argtup, star_indx, rv, ew, λ0, sigma, injectindx, injectfiber = ovtup
        ival = argtup[1]
        intup = argtup[2:end]
        (release_dir,redux_ver,tele,field,plate,mjd,fiberindx) = intup

        teleind = (tele[1:6] == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiberindx
        # Get Throughput Fluxing 
        fluxingcache = cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir)
        if !isfile(fluxingcache)
            dirName = splitdir(fluxingcache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            getAndWrite_fluxing(release_dir,redux_ver,tele,field,plate,mjd,cache_dir=cache_dir)
        end
        
        starname_old = cache_starname(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
        if (isfile(starname_old) & caching)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(starname_old)
            starscale,framecnts,a_relFlux,b_relFlux,c_relFlux,cartVisit,ingest_bit = metaexport
        else
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
            starscale,framecnts,a_relFlux,b_relFlux,c_relFlux,cartVisit,ingest_bit = metaexport
            if caching
                dirName = splitdir(starname_old)[1]
                if !ispath(dirName)
                    mkpath(dirName)
                end
                serialize(starname_old,[fvec, fvarvec, cntvec, chipmidtimes, metaexport])
            end
        end
        
        outup = (release_dir,redux_ver,tele*lpad(injectindx,7,"0")*"i",field,plate,mjd,injectfiber)
        starcache = cache_starname(tele*lpad(injectindx,7,"0")*"i",field,plate,mjd,injectfiber,cache_dir=cache_dir,inject_cache_dir=inject_cache_dir)
        
        if !isfile(starcache)
            starcomp = x_starContinuum[:,star_indx].*(1 .+ x_starLines[:,star_indx])

            # analytic shift, so don't shift the LSF
            Ksp = if adjfibindx>300
                deserialize(prior_dict["LSF_mat_LCO"]*"6_"*lpad(adjfibindx-300,3,"0")*".jdat");
            else
                deserialize(prior_dict["LSF_mat_APO"]*"6_"*lpad(adjfibindx,3,"0")*".jdat");
            end
            nvecLSF = dropdims(sum(Ksp,dims=2),dims=2);
            for (dib_ind, dib_center_lambda) in enumerate(dib_center_lambda_lst)
                dibcomp = Ksp*gauss1d_ew(ew[dib_ind],λ0[dib_ind],sigma[dib_ind],x_model)./nvecLSF
                starcomp .*= (1 .+ dibcomp)
            end
            starcomp .-= 1 # this is now the negative flux we should add the observation to make the simulated observation

            outfvec = fvec.+starcomp
            rng = MersenneTwister(mjd+adjfibindx) # seed on mjd and fiber index so unique

            outcntvec = cntvec

            framecnts = maximum(outcntvec)
            simplemsk = (outcntvec.==framecnts)
            
            starscale_o = if count(simplemsk .& (.!isnan.(outfvec)))==0
                NaN
            else
                abs(nanzeromedian(outfvec[simplemsk]))
            end

            outvar = fvarvec # assume DIB does not impact noise model significantly

            # just passing along the meta data
            metaexport = (starscale_o,framecnts,a_relFlux,b_relFlux,c_relFlux,cartVisit,ingest_bit)

            dirName = splitdir(starcache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(starcache,[outfvec, outvar, outcntvec, chipmidtimes, metaexport])
            return (injectindx, outup...), starcache
        else
            return (injectindx, outup...), starcache
        end
    end
end

## Generate Injection Parameters
rng = MersenneTwister(rnd_seed)
injectfiber = fiberindx*ones(Int,nsamp) # this should not be adjusted
injectindx = 1:nsamp
star_indx = rand(rng,obs_indices2use,nsamp);
star_tup = ntuplst[star_indx]
dib_sig = []
dib_lam = []
dib_ew = []
for dib_center_lambda in dib_center_lambda_lst
    push!(dib_sig,rand(rng,Uniform(dib_sig_range...),nsamp));
    push!(dib_lam,rand(rng,Uniform(dib_center_lambda*(1+dib_vel_range[1]/c),dib_center_lambda*(1+dib_vel_range[2]/c)),nsamp));
    if dibs_on
        push!(dib_ew,rand(rng,Uniform(dib_ew_range...),nsamp))
    else
        push!(dib_ew,zeros(nsamp))
    end
end

itarg = Iterators.zip(star_tup,star_indx,eachrow(hcat(dib_ew...)),eachrow(hcat(dib_lam...)),eachrow(hcat(dib_sig...)),injectindx,injectfiber);

## Save Injection Parameters to Disk
fname = prior_dict["out_dir"]*"inject_params.h5"
dirName = splitdir(fname)[1]
if !ispath(dirName)
    mkpath(dirName)
end
f = h5open(fname,"w")
write(f,"injectfiber",injectfiber)
write(f,"injectindx",collect(injectindx))
write(f,"star_indx",star_indx)
for (dibind, dib_center_lambda) in enumerate(dib_center_lambda_lst)
    write(f,"dib_sig_$(dib_center_lambda)",dib_sig[dibind])
    write(f,"dib_lam_$(dib_center_lambda)",dib_lam[dibind])
    write(f,"dib_ew_$(dib_center_lambda)",dib_ew[dibind])
end
close(f)

## Create an injection test series
@everywhere take_draw_partial(ovtup) = take_draw(ovtup,dib_center_lambda_lst=dib_center_lambda_lst,inject_cache_dir=prior_dict["inject_cache_dir"],cache_dir=prior_dict["local_cache"])
pout = @showprogress pmap(take_draw_partial,itarg);

## Write out runlist needed to run the apMADGICS.jl on it
input_lst = []
starcache_lst = []
for (subindx, subout) in enumerate(pout)
    push!(input_lst,subout[1])
    push!(starcache_lst,subout[2])
end

println("All of the injected spectra wrote to unique cache files: ", length(starcache_lst) == length(unique(starcache_lst)))

serialize(prior_dict["out_dir"]*"injection_input_lst_"*lpad(adjfibindx,3,"0")*".jdat",input_lst);
