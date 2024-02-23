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

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    runlist_range = 335 #1:600 #295, 245, 335, 101

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Input List (not really a prior, but an input file we search for stars conditioned on)
    prior_dict["tellStdGood"] = prior_dir*"2024_02_21/outlists/summary/dr17_dr17_good_tell_list.txt"
    prior_dict["tell_runlist"] = prior_dir*"2024_02_21/outlists/tell/dr17_dr17_tell_input_lst_plate_msked_"

    # Data for Detector Cals (not really a prior, but an input the results depend on in detail)
    prior_dict["chip_fluxdep_err_correction"] = src_dir*"data/chip_fluxdep_err_correction.jdat"
    prior_dict["medframes_APO"] = src_dir*"data/medframes_apo.jdat" # last made 2023_04_01 by AKS
    prior_dict["medframes_LCO"] = src_dir*"data/medframes_lco.jdat" # last made 2023_04_07 by AKS
    prior_dict["LSF_mat_APO"] = prior_dir0*"2023_04_01/mat_lsf_out/sp_combolsfmat_norm_6_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir0*"2023_04_07/mat_lsf_out/sp_combolsfmat_norm_6_" # last made 2023_04_07 by AKS

    # Sky tellDiv priors
    sky_base = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/"
    prior_dict["skycont_tellDiv"] = sky_base*"sky_priors/APOGEE_skycont_tellDiv_svd_30_f"
    prior_dict["skyline_tellDiv"] = sky_base*"sky_priors/APOGEE_skyline_faint_tellDiv_GSPICE_svd_120_f"
end

@everywhere begin
    # I should revisit the error bars in the context of chi2 versus frame number trends
    global err_correct_Dict = deserialize(prior_dict["chip_fluxdep_err_correction"])

    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg)
    x_model = 15000:0.01:17000

    c = 299792.458; # in km/s
    delLog = 6e-6; 

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
    telvec = zeros(length(wavetarg));
end

@everywhere begin
    bb10_000 = blackbody.(Ref(10_000), x_model*1e-8);

    brackett = [
        16811.129, #11
        16411.692, #12 
        16113.731, #13
        15884.898, #14
        15704.969, #15
        15560.716, #16
        15443.156, #17
        15345.999, #18
        # 15264.725, #19
        # 15196.013, #20
    ]

    amplst = [
        -5, #11
        -4, #12 
        -4.5, #13
        -2, #14
        -2, #15
        -1.5, #16
        -1.5, #17
        -0.2, #18
        -0.2, #19
        -0.2, #20
    ]

    nseries = 11:20;

    bsigrng = 10:1:50
    sigrng = bsigrng
    sigstep0 = step(sigrng)
    minoffset0  = round(Int,sigrng[1] / sigstep0)-1

    Vbrack = zeros(length(wavetarg),length(brackett),length(bsigrng))
    for (cind,cent) in enumerate(brackett)
        for (sigind,sig) in enumerate(bsigrng)
            Vbrack[:,cind,sigind].=lorentzian.(wavetarg,Ref(2*amplst[cind]),Ref(cent),Ref(sig))
        end
    end
end

@everywhere begin
    function sky_res_extract(outvec,outvar,simplemsk,V_skyline_faint,V_skycont)
        ## Select data for use (might want to handle mean more generally)
        wave_obs = wavetarg[simplemsk]
        Xd_obs = outvec[simplemsk];

        ## Set up residuals prior
        A = Diagonal(outvar[simplemsk]);
        Ainv = Diagonal(1 ./outvar[simplemsk]);

        ## Set up priors
        V_skyline_c = V_skyline_faint
        V_skyline_r = V_skyline_c[simplemsk,:]

        V_skycont_c = V_skycont
        V_skycont_r = V_skycont_c[simplemsk,:]

        # Compute sky line/continuum separation
        Vcomb = hcat(V_skyline_r,V_skycont_r)
        Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
        x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, 
            (A, V_skycont_r, ),
            (A, V_skycont_c, ))
        return x_comp_lst
    end

    # function genModSamp(intup)
    #     Teff,Av,Rv,Tfunindx,Tfracindx = intup
    #     bbs = blackbody.(Ref(Teff), x_model*1e-8);
    #     rvec = redden_mult(x_model,Av,Rv);
    #     return tellFracSamples[:,Tfracindx].*TfunSamples[:,Tfunindx].*((Ksp*(rvec.*bbs))./nvecLSF);
    # end

    sigScanFun(x; step = sigstep0, minoffset = minoffset0)  = round(Int,x / step) .-minoffset

    function chi2_wrapper(sval,intup)
        (simplemsk,Ctotinv,Xd_obs,Dscale,V_new,pre_Vslice) = intup
        # transform val to index
        rval = indInt(sval)
        tval = indTenth(sval)
        pre_Vslice .= view(ShiftedArrays.circshift(V_new,(rval,0)),simplemsk,:)
        pre_Vslice .*= Dscale 
        return woodbury_update_inv_tst(
            Ctotinv,
            Xd_obs,
            pre_Vslice
        )
    end

    function chi2_wrapper_b(sval,intup)
        (simplemsk,Ctotinv,Xd_obs,Dscale,V_new,pre_Vslice,rval) = intup
        sigindx = sigScanFun(sval)
        pre_Vslice .= view(ShiftedArrays.circshift(view(V_new,:,sigindx),(rval,0)),simplemsk,:)
        pre_Vslice .*= Dscale 
        return woodbury_update_inv_tst(
            Ctotinv,
            Xd_obs,
            pre_Vslice
        )
    end

    function cheb_fit_refine(outvec,outvar,simplemsk,Vpoly_scaled,Vbrack)
        ## Select data for use (might want to handle mean more generally)
        wave_obs = wavetarg[simplemsk]
        Xd_obs = outvec[simplemsk];

        ## Set up residuals prior
        A = Diagonal(outvar[simplemsk]);
        Ainv = Diagonal(1 ./outvar[simplemsk]);

        ## Set up priors
        V_smooth_c = Vpoly_scaled
        V_smooth_r = V_smooth_c[simplemsk,:]

        Vcomb_cur = V_smooth_r
        Ctotinv_cur = LowRankMultMat([Ainv,Vcomb_cur],wood_precomp_mult,wood_fxn_mult);
        x_comp_lst = deblend_components_all_asym(Ctotinv_cur, Xd_obs, 
            (V_smooth_r, ), 
            (V_smooth_r, ))

        Vbrack_red = Vbrack[:,:,21]
        pre_Vslice = zeros(count(simplemsk),size(Vbrack_red,2));
        chi2_wrapper_partial = Base.Fix2(chi2_wrapper,(simplemsk,Ctotinv_cur,Xd_obs,x_comp_lst[1],Vbrack_red,pre_Vslice))

        lvl1 = -70:1//1:70
        lout = sampler_1d_dense_var(chi2_wrapper_partial,lvl1,stepx=1)
        inshift = lout[1][3]

        lvl1d = bsigrng
        pre_Vslice = zeros(count(simplemsk),1);
        V_brackett_all_r = zeros(count(simplemsk),length(brackett))
        V_brackett_all_c = zeros(length(wavetarg),length(brackett))
        for (cind,cent) in enumerate(brackett)
            Vin = ShiftedArrays.circshift(view(Vbrack,:,cind,:),(inshift,0))
            chi2_wrapper_partial = Base.Fix2(chi2_wrapper_b,(simplemsk,Ctotinv_cur,Xd_obs,x_comp_lst[1],Vin,pre_Vslice,inshift))
            lout = sampler_1d_dense_var(chi2_wrapper_partial,lvl1d,stepx=1)

            V_brackett_all_c[:,cind] = Vin[:,lout[1][5]]
            V_brackett_all_r[:,cind] = V_brackett_all_c[simplemsk,cind]
            V_brackett_all_r[:,cind] .*= x_comp_lst[1]

            Vcomb_fut = hcat(Vcomb_cur,V_brackett_all_r[:,cind]);
            Ctotinv_fut = LowRankMultMat([Ainv,Vcomb_fut],wood_precomp_mult,wood_fxn_mult);
            Ctotinv_cur, Ctotinv_fut = Ctotinv_fut, Ctotinv_cur; Vcomb_cur, Vcomb_fut = Vcomb_fut, Vcomb_cur
        end

        x_comp_lst = deblend_components_all_asym(Ctotinv_cur, Xd_obs, 
            (A, V_smooth_r), 
            (A, V_smooth_c))

        return x_comp_lst[2]
    end

    function ingest_tellVisit_stack(argtup;cache_dir="./local_cache_Tfun")
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        # Get Throughput Fluxing 
        fluxingcache = cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir)
        if !isfile(fluxingcache)
            dirName = splitdir(fluxingcache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            getAndWrite_fluxing(release_dir,redux_ver,tele,field,plate,mjd,cache_dir=cache_dir)
        end

        skyLineCache_tellDiv = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
        if !isfile(skyLineCache_tellDiv)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvec = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
            dirName = splitdir(skyLineCache_tellDiv)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(skyLineCache_tellDiv,[fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvec])
        end
    end

    function extract_Tfun_wrapper(argtup,kout,Vpoly_scaled,Vbrack,V_skyline_faint,V_skycont,skymsk;cache_dir="./local_cache_Tfun")
        ival = argtup[1]
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache_tellDiv = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
        fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvecN = deserialize(skyLineCache_tellDiv)
        simplemsk = (cntvec.==maximum(cntvec)) .& skymsk;

        res_r, skycont_c = sky_res_extract(fvec,fvarvec,simplemsk,V_skyline_faint,V_skycont)

        datvec = nanify(res_r,simplemsk) .+ skycont_c
        BBscaled = (nanmean(datvec)/mean(kout)).*kout;
        testvec = datvec./BBscaled;
        testvar = fvarvec./(BBscaled.^2)

        tfun_c = cheb_fit_refine(testvec,testvar,simplemsk,Vpoly_scaled,Vbrack)

        return tfun_c
    end
end

@everywhere begin
    function get_Tfun_samples(adjfibindx;loc_parallel=false)

        savename = "tell_prior_disk/tfun_samples_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            ntuplst = deserialize(prior_dict["tell_runlist"]*lpad(adjfibindx,3,"0")*".jdat")

            medframes = if adjfibindx>300
                deserialize(prior_dict["medframes_LCO"]);
            else
                deserialize(prior_dict["medframes_APO"]);
            end

            Ksp = if adjfibindx>300
                deserialize(prior_dict["LSF_mat_LCO"]*lpad(adjfibindx-300,3,"0")*".jdat");
            else
                deserialize(prior_dict["LSF_mat_APO"]*lpad(adjfibindx,3,"0")*".jdat");
            end

            #global nvecLSF = dropdims(sum(Ksp,dims=2),dims=2); # used only in starCont sample gen
            kout = Ksp*bb10_000;

            f = h5open(prior_dict["skycont_tellDiv"]*lpad(adjfibindx,3,"0")*".h5")
            V_skycont = read(f["Vmat"])
            chebmsk_exp = convert.(Bool,read(f["chebmsk_exp"]))
            close(f)
        
            f = h5open(prior_dict["skyline_tellDiv"]*lpad(adjfibindx,3,"0")*".h5")
            V_skyline_faint = read(f["Vmat"])
            submsk_faint = convert.(Bool,read(f["submsk"]))
            close(f)
            
            skymsk = chebmsk_exp .& submsk_faint;

            Vpoly_scaled, msknall = generate_poly_prior(adjfibindx,medframes)

            extract_Tfun_wrapper_bound(ival) = extract_Tfun_wrapper(ival,kout,Vpoly_scaled,Vbrack,V_skyline_faint,V_skycont,skymsk)

            # Write Tell Obs to Disk
            pout = if loc_parallel
                @showprogress pmap(ingest_tellVisit_stack,ntuplst);
            else
                map(ingest_tellVisit_stack,ntuplst);
            end

            pout = if loc_parallel
                @showprogress pmap(extract_Tfun_wrapper_bound,ntuplst);
            else
                map(extract_Tfun_wrapper_bound,ntuplst);
            end

            TfunSamples = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                TfunSamples[:,i].=pout[i]
            end

            dirName = splitdir(savename)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(savename,TfunSamples)
        end
    end
end

get_Tfun_samples(runlist_range)
# good_tell_list = convert.(Int,readdlm("../2024_02_21/outlists/summary/dr17_dr17_good_tell_list.txt",','))[:,1]
# @showprogress pmap(get_Tfun_samples,good_tell_list)