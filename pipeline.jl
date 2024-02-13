## This is the main pipeline that will batch over APOGEE files
# Author - Andrew Saydjari, CfA

import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("./"); Pkg.instantiate(); Pkg.precompile()
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
using MKL
using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager(),exeflags=["--project=./"])
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker allocation took $dt"); t_then = t_now; flush(stdout)

@everywhere begin
    using MKL
    using LinearAlgebra
    BLAS.set_num_threads(1)
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using ThreadPinning

    prior_dir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/"
    src_dir = "./"
    include(src_dir*"src/utils.jl")
    include(src_dir*"src/gridSearch.jl")
    include(src_dir*"src/componentAndPosteriors.jl")
    include(src_dir*"src/fileNameHandling.jl")
    include(src_dir*"src/ingest.jl")
    include(src_dir*"src/lowRankPrescription.jl")
    include(src_dir*"src/marginalizeEW.jl")
    include(src_dir*"src/spectraInterpolation.jl")
    include(src_dir*"src/chi2Wrappers.jl")
    
    using StatsBase, ProgressMeter
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

# These global allocations for the injest are messy... but we plan on changing the ingest
# relatively soon... so don't worry for now.
# SHOULD probably make a flag for using DD correction or not (do once we settle on the formalism)
@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg)
    
    c = 299792.458; # in km/s
    delLog = 6e-6; 
    # pixscale = (10^(delLog)-1)*c;

    # nothing to do on size here, if anything expand
    f = h5open(prior_dir*"2023_07_22/dib_priors/precomp_dust_1_analyticDeriv_stiff.h5")
    global V_dib_noLSF = read(f["Vmat"])
    close(f)

    f = h5open(prior_dir*"2023_07_22/dib_priors/precomp_dust_3_analyticDeriv_soft.h5")
    global V_dib_noLSF_soft = read(f["Vmat"])
    close(f)

    alpha = 1;
    f = h5open(prior_dir*"2023_08_22/starLine_priors/APOGEE_stellar_kry_50_subpix_th22500.h5")
    global V_subpix_refLSF = alpha*read(f["Vmat"])
    close(f)

    Xd_stack = zeros(3*2048)
    Xd_std_stack = zeros(3*2048)
    waveobs_stack = zeros(3*2048)
    waveobs_stack_old = zeros(3*2048)
    pixmsk_stack = zeros(Int,3*2048)
    fullBit = zeros(Int,3*2048);
    outvec = zeros(length(wavetarg))
    outvar = zeros(length(wavetarg))
    cntvec = zeros(Int,length(wavetarg));
end

# it would be great to move this into a parameter file that is read for each run
@everywhere begin
    refine_iters = 1

    # Moon Wave
    mlvl1 = -2:1//10:2
    mlvl_tuple = (mlvl1,)
    # tuple1dprint(mlvl_tuple)

    # Star Wave
    lvl1 = -70:1//2:70
    lvl2 = -8:2//10:8
    lvl3 = -3:1//10:3
    slvl_tuple = (lvl1,lvl2,lvl3)
    # tuple1dprint(slvl_tuple)

    # (Wave, Sig) DIB
    dib_center_lst = [15273, 15273]#, 15653]
    lvl1d = ((-150:4:150),(18//10:18//10))
    lvl2d = ((0:0), (-7//5:4//100:11//5))
    lvl3d = ((-18:2//10:18), (0:0))
    lvl4d = ((0:0), (-90//100:2//100:90//100))
    lvl5d = ((-1:2//10:1), (0:0))
    lvl6d = ((0:0), (-10//100:2//100:10//100))
    lvl7d = ((-6//10:2//10:6//10), (0:0))
    lvl8d = ((0:0), (-6//100:2//100:6//100))
    lvl9d = ((-4//10:1//10:4//10), (-4//100:1//100:4//100));
    lvltuple = (lvl1d, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    lvl1d_2 = ((-60:4:60),(18//10:18//10))
    lvltuple_2 = (lvl1d_2, lvl2d, lvl3d, lvl4d, lvl5d, lvl6d, lvl7d, lvl8d, lvl9d);
    lvltuple_lst = [lvltuple, lvltuple_2]
    # tuple2dprint(lvltuple)

    # Flux marginalize region
    sigMarg0 = -50//100:10//100:50//100
    svalMarg0 = -0//10:1//10:0//10;
end

@everywhere begin
    function pipeline_single_spectra(argtup; caching=true, sky_caching=true, skyCont_off=false, skyLines_off=false, rv_chi2res=false, rv_split=true, cache_dir="../local_cache", inject_cache_dir=prior_dir*"2024_02_08/inject_local_cache")
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]
        out = []

        # Get Throughput Fluxing 
        fluxingcache = cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir)
        if !isfile(fluxingcache)
            dirName = splitdir(fluxingcache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            getAndWrite_fluxing(release_dir,redux_ver,tele,field,plate,mjd,cache_dir=cache_dir)
        end

        # Get Sky Prior
        skycache = cache_skyname(tele,field,plate,mjd,cache_dir=cache_dir)
        if (isfile(skycache) & sky_caching)
            meanLocSky, VLocSky = deserialize(skycache)
        else
            meanLocSky, VLocSky = getSky4visit(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,skymsk,V_skyline_bright,V_skyline_faint,V_skycont,caching=sky_caching,cache_dir=cache_dir)
            if sky_caching
                dirName = splitdir(skycache)[1]
                if !ispath(dirName)
                    mkpath(dirName)
                end
                serialize(skycache,[meanLocSky, VLocSky])
            end
        end

        # Get the Star
        starcache = cache_starname(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir,inject_cache_dir=inject_cache_dir)
        if (isfile(starcache) & caching)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(starcache)
            starscale,framecnts,varoffset,varflux,a_relFlux,b_relFlux,c_relFlux,cartVisit = metaexport
        elseif tele[end]=='i'
            @warn "Injections not found at injection cache dir!"
        else
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
            starscale,framecnts,varoffset,varflux,a_relFlux,b_relFlux,c_relFlux,cartVisit = metaexport
            if caching
                dirName = splitdir(starcache)[1]
                if !ispath(dirName)
                    mkpath(dirName)
                end
                serialize(starcache,[fvec, fvarvec, cntvec, chipmidtimes, metaexport])
            end
        end
        simplemsk = (cntvec.==framecnts) .& skymsk;
        
        push!(out,(count(simplemsk), starscale, framecnts, chipmidtimes, varoffset, varflux, a_relFlux, b_relFlux, c_relFlux, cartVisit, nanify(fvec[simplemsk],simplemsk), nanify(fvarvec[simplemsk],simplemsk))) # 1

        if skyCont_off
            meanLocSky.=0
            VLocSky.=0
        end
        if skyLines_off
            V_skyline_bright.=0
            V_skyline_faint.=0
        end
        ## Select data for use (might want to handle mean more generally)
        Xd_obs = (fvec.-meanLocSky)[simplemsk]; #I think an outvec to fvec here was the key caching issue
        wave_obs = wavetarg[simplemsk]

        ## Set up residuals prior
        A = Diagonal(fvarvec[simplemsk]);
        Ainv = Diagonal(1 ./fvarvec[simplemsk]);

        ## Create interpolation matrices
        fullBitprox = zeros(Int,length(wavetarg))
        fullBitprox[.!simplemsk] .= 2^4

        ## Set up priors
        # V_skyline_bright_c = V_skyline_bright
        # V_skyline_bright_r = V_skyline_bright_c[simplemsk,:]
        V_skyline_faint_c = V_skyline_faint
        V_skyline_faint_r = V_skyline_faint_c[simplemsk,:]
        # V_skyline_tot_c = hcat(V_skyline_bright_c,V_skyline_faint_c)
        # V_skyline_tot_r = hcat(V_skyline_bright_r,V_skyline_faint_r)
        V_skyline_tot_c = V_skyline_faint_c
        V_skyline_tot_r = V_skyline_faint_r
        V_locSky_c = VLocSky
        V_locSky_r = V_locSky_c[simplemsk,:]
        V_starCont_c = starscale*V_starcont
        V_starCont_r = V_starCont_c[simplemsk,:]

        ## Solve RV of Star
        # compute stellar continuum to modify stellar line prior
        Vcomb_skylines = hcat(V_skyline_tot_r,V_locSky_r,V_starCont_r);
        Ctotinv_skylines = LowRankMultMatIP([Ainv,Vcomb_skylines],wood_precomp_mult_mat([Ainv,Vcomb_skylines],(size(Ainv,1),size(V_subpix,2))),wood_fxn_mult,wood_fxn_mult_mat!);
        x_comp_lst = deblend_components_all(Ctotinv_skylines, Xd_obs, (V_starCont_r,))
        starCont_Mscale = x_comp_lst[1]

        # now take out the skylines to be included in the scanning
        Vcomb_cur = hcat(V_locSky_r,V_starCont_r);
        Ctotinv_cur = LowRankMultMatIP([Ainv,Vcomb_cur],wood_precomp_mult_mat([Ainv,Vcomb_cur],(size(Ainv,1),size(V_subpix,2))),wood_fxn_mult,wood_fxn_mult_mat!);

        # compute delta chi2 for adding skylines (helps normalize the joint chi2 below with starLines)
        chi2skyoffset = woodbury_update_inv_tst(
            LowRankMultMatIP([Ainv,Vcomb_cur],wood_precomp_mult_mat([Ainv,Vcomb_cur],(size(Ainv,1),size(V_skyline_tot_r,2))),wood_fxn_mult,wood_fxn_mult_mat!),
            Xd_obs,
            V_skyline_tot_r
        )
  
        pre_Vslice = zeros(count(simplemsk),size(V_subpix,2))
        chi2_wrapper_partial = if rv_chi2res
            Base.Fix2(chi2_wrapper_res,(simplemsk,Ctotinv_cur,Xd_obs,starCont_Mscale,V_subpix,pre_Vslice,A))
        elseif rv_split
            AinvV1 = Ctotinv_cur*V_skyline_tot_r
            XdAinvV1 = reshape(Xd_obs,1,:)*AinvV1
            V1TAinvV1 = V_skyline_tot_r'*AinvV1
            Base.Fix2(chi2_wrapper_split,(simplemsk,Ctotinv_cur,Xd_obs,starCont_Mscale,V_subpix,pre_Vslice,AinvV1,XdAinvV1,V1TAinvV1,chi2skyoffset))
        else
            Base.Fix2(chi2_wrapper,(simplemsk,Ctotinv_cur,Xd_obs,starCont_Mscale,V_subpix,pre_Vslice))
        end
        lout = sampler_1d_hierarchy_var(chi2_wrapper_partial,slvl_tuple,minres=1//10,stepx=8)
        push!(out,lout) # 2

        # update the Ctotinv to include the stellar line component (iterate to refine starCont_Mscale)
        svalc = lout[1][3]
        for i=1:refine_iters
            Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],simplemsk,starCont_Mscale,Vcomb_skylines,V_subpix,V_subpix_refLSF)
            x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,))
            starCont_Mscale = x_comp_lst[1]
        end
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],simplemsk,starCont_Mscale,Vcomb_skylines,V_subpix,V_subpix_refLSF)
        
        # do a component save without the 15273 DIB
        x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
            (A, V_skyline_faint_r, V_locSky_r, V_starCont_r, V_starlines_r, V_starlines_r),
            (A, V_skyline_faint_r, V_locSky_r, V_starCont_r, V_starlines_c, I),
        )
        push!(out,x_comp_lst[1]'*(Ainv*x_comp_lst[1])) # 3
        x_comp_out = [nanify(x_comp_lst[1],simplemsk)./sqrt.(fvarvec), nanify(x_comp_lst[1],simplemsk), 
                        # nanify(x_comp_lst[2][skymsk_bright[simplemsk]],simplemsk .& skymsk_bright), nanify(x_comp_lst[3][skymsk_faint[simplemsk]],simplemsk .& skymsk_faint), 
                        nanify(x_comp_lst[2][skymsk_faint[simplemsk]],simplemsk .& skymsk_faint), 
                        nanify(x_comp_lst[3].+meanLocSky[simplemsk],simplemsk), nanify(x_comp_lst[4],simplemsk),
                        x_comp_lst[5:end]..., nanify((fvec[simplemsk].-(x_comp_lst[2].+x_comp_lst[3].+meanLocSky[simplemsk]))./ x_comp_lst[4],simplemsk)]
        push!(out,x_comp_out) # 4
        dflux_starlines = sqrt_nan.(get_diag_posterior_from_prior_asym(Ctotinv_fut, V_starlines_c, V_starlines_r))
        push!(out,dflux_starlines) # 5
                
        # prepare multiplicative factors for DIB prior
        x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,V_starlines_r))
        starCont_Mscale = x_comp_lst[1]
        starFull_Mscale = x_comp_lst[1].+x_comp_lst[2]
        
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_skylines.matList[1],simplemsk,starCont_Mscale,Vcomb_skylines,V_subpix,V_subpix_refLSF)
        Ctotinv_cur, Ctotinv_fut = Ctotinv_fut, Ctotinv_cur; Vcomb_cur, Vcomb_fut = Vcomb_fut, Vcomb_cur # swap to updated covariance finally
        
        # currently, this is modeling each DIB seperately... I think we want to change this later, just easier parallel structure
        pre_Vslice = zeros(count(simplemsk),size(V_dib,2))
        for dib_ind = 1:length(dib_center_lst) # eventually need to decide if these are cumulative or not
            lvltuple_dib = lvltuple_lst[dib_ind]
            dib_center = dib_center_lst[dib_ind]
            scan_offset = findmin(abs.(wavetarg.-dib_center_lst[dib_ind]))[2].-findmin(abs.(wavetarg.-dib_center_lst[1]))[2]
            
            ## Solve DIB parameters (for just 15273), not any more, just a single DIB
            # one of the main questions is how many time to compute components and where
            chi2_wrapper_partial = Base.Fix2(chi2_wrapper2d,(simplemsk,Ctotinv_cur,Xd_obs,wave_obs,starFull_Mscale,Vcomb_cur,V_dib,pre_Vslice,dib_center,scan_offset))
            lout = sampler_2d_hierarchy_var(chi2_wrapper_partial,lvltuple_dib,step1=3,step2=3,minres1=1//10,minres2=1//100)
            opt_tup = lout[1][3]
            push!(out,lout) # 6

            ## Shift the marginalization sampling (should this be wrapped inside the function?)
            # especially because we need to do bounds handling
            svalMarg = svalMarg0 .+ opt_tup[1]
            sigMarg = shift_trim_range(sigMarg0,opt_tup[2]; minv=4//10, maxv=4)
            samp_lst = Iterators.product(svalMarg,sigMarg)

            intupf = (simplemsk,Ctotinv_cur,Xd_obs,wave_obs,starFull_Mscale,Vcomb_cur,V_dib,pre_Vslice,dib_center,scan_offset)
            chi2lst, fluxlst, dfluxlst = sample_chi2_flux_dflux(samp_lst,intupf) #shouldn't this take chi2_wrapper_partial as an argument?
            refchi2val = minimum(chi2lst) #this should just be set to the min found at the 2d step
            lout = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
            push!(out,lout) # 7

            # Compute some final components for export (still need to implement DIB iterative refinement)
            Ctotinv_fut, Vcomb_fut, V_dibc, V_dibr = update_Ctotinv_Vdib_asym(
                opt_tup,Ctotinv_cur.matList[1],simplemsk,starFull_Mscale,Vcomb_cur,V_dib_soft,V_dib_noLSF_soft,scan_offset)

            x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
                (A, V_skyline_faint_r, V_locSky_r, V_starCont_r, V_starlines_r, V_dibr),
                (A, V_skyline_faint_r, V_locSky_r, V_starCont_r, V_starlines_c, V_dibc),
            )
            push!(out,x_comp_lst[1]'*(Ainv*x_comp_lst[1])) # 8
            # I am not sure that during production we really want to run and output full sets of components per DIB
            # I would like to fill NaNs in chip gaps for the sky/continuum components
            # revisit that when we revisit the interpolations before making other fiber priors
            x_comp_out = [nanify(x_comp_lst[1],simplemsk)./sqrt.(fvarvec), nanify(x_comp_lst[1],simplemsk),
                        # nanify(x_comp_lst[2][skymsk_bright[simplemsk]],simplemsk .& skymsk_bright), nanify(x_comp_lst[3][skymsk_faint[simplemsk]],simplemsk .& skymsk_faint), 
                        nanify(x_comp_lst[2][skymsk_faint[simplemsk]],simplemsk .& skymsk_faint), 
                        nanify(x_comp_lst[3].+meanLocSky[simplemsk],simplemsk), nanify(x_comp_lst[4],simplemsk),
                        x_comp_lst[5:end]...]

            push!(out,x_comp_out) # 9
        end
                        
        return out
    end
end

@everywhere begin
    function multi_spectra_batch(indsubset; out_dir="../outdir", ddstaronly=false)
        ### Set up
        out = []
        startind = indsubset[1][1]
        global tele = indsubset[1][4]
        fiberindx = indsubset[1][end]
        teleind = (tele[1:6] == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiberindx

        ### Save and cache restart handling
        savename = join([out_dir,lpad(adjfibindx,3,"0"),"apMADGICS_fiber_"*lpad(adjfibindx,3,"0")*"_batch_"*lpad(startind,7,"0")*".h5"],"/")
        dirName = splitdir(savename)[1]
        if !ispath(dirName)
            mkpath(dirName)
        end
        if !isfile(savename)
            prior_load_needed = if @isdefined loaded_adjfibindx
                if adjfibindx != loaded_adjfibindx
                    true
                else
                    false
                end
            else
                true
            end
            if prior_load_needed
                ### Need to load the priors here
                f = h5open(prior_dir*"2024_01_31/sky_priors/APOGEE_skycont_svd_30_f"*lpad(adjfibindx,3,"0")*".h5")
                global V_skycont = read(f["Vmat"])
                chebmsk_exp = convert.(Bool,read(f["chebmsk_exp"]))
                close(f)

                f = h5open(prior_dir*"2024_02_04/sky_priors/APOGEE_skyline_bright_svd_120_20_8_6_1en3_f"*lpad(adjfibindx,3,"0")*".h5") #revert temp
                global V_skyline_bright = read(f["Vmat"])
                submsk_bright = convert.(Bool,read(f["submsk"]))
                close(f)

                f = h5open(prior_dir*"2024_02_04/sky_priors/APOGEE_skyline_faint_svd_120_20_8_6_1en3_f"*lpad(adjfibindx,3,"0")*".h5") #revert temp
                global V_skyline_faint = read(f["Vmat"])
                submsk_faint = convert.(Bool,read(f["submsk"]))
                close(f)

                global skymsk_bright = chebmsk_exp .& submsk_bright #.& msk_starCor;
                global skymsk_faint = chebmsk_exp .& submsk_faint #.& msk_starCor;
                # global skymsk = chebmsk_exp .& (submsk_bright .| submsk_faint) #.& msk_starCor;
                global skymsk = chebmsk_exp .& submsk_faint # completely masking all bright lines;

                f = h5open(prior_dir*"2023_07_22/star_priors/APOGEE_starcont_svd_60_f"*lpad(adjfibindx,3,"0")*".h5")
                global V_starcont = read(f["Vmat"])
                close(f)

                # can consider changing dimension at the full reduction stage
                f = h5open(prior_dir*"2023_08_22/starLine_priors/APOGEE_stellar_kry_50_subpix_"*lpad(adjfibindx,3,"0")*".h5")
                global V_subpix = alpha*read(f["Vmat"])
                close(f)
                if ddstaronly
                    global V_subpix_refLSF = V_subpix
                end

                f = h5open(prior_dir*"2023_07_22/dib_priors/precomp_dust_1_analyticDerivLSF_stiff_"*lpad(adjfibindx,3,"0")*".h5")
                global V_dib = read(f["Vmat"])
                close(f)

                f = h5open(prior_dir*"2023_07_22/dib_priors/precomp_dust_3_analyticDerivLSF_soft_"*lpad(adjfibindx,3,"0")*".h5")
                global V_dib_soft = read(f["Vmat"])
                close(f)
                GC.gc()
            end
            global loaded_adjfibindx = adjfibindx
            
            ### Single spectrum loop
            for (ind,indval) in enumerate(indsubset)
                push!(out,pipeline_single_spectra(indval; caching=true))
            end

            ### Save Exporting
            metai = 1
            RVind, RVchi, RVcom, strpo = 2, 3, 4, 5
            DIBind, EWind, DIBchi, DIBcom = 6, 7, 8, 9
            dibsavesz = 4

            ## RV Block
            RVextract = [
                # meta info
                (x->x[metai][1],                        "data_pix_cnt"),
                (x->x[metai][2],                        "starscale"),
                (x->x[metai][3],                        "frame_counts"),
                (x->x[metai][4],                        "chip_midtimes"),
                (x->x[metai][5],                        "varoffset"),
                (x->x[metai][6],                        "varflux"),
                (x->x[metai][7],                        "a_relFlux"),
                (x->x[metai][8],                        "b_relFlux"),
                (x->x[metai][9],                        "c_relFlux"),
                # (x->x[metai][10],                       "cartVisit"),
                (x->x[metai][11],                       "flux"),
                (x->x[metai][12],                       "fluxerr2"),
                (x->adjfibindx,                         "adjfiberindx"),

                (x->Float64.(x[RVind][1][1]),           "RV_pixoff_final"),
                (x->Float64.(x[RVind][1][3]),           "RV_pixoff_disc_final"),
                (x->x[RVind][1][2],                     "RV_minchi2_final"),
                (x->x[RVind][1][6],                     "RV_flag"),
                (x->x[RVind][1][7],                     "RV_pix_var"),
                                    
                (x->x[RVchi][1],                        "RVchi2_residuals"),
                                    
                (x->x[RVind][2][1][3],                  "RV_p5delchi2_lvl1"),
                (x->x[RVind][2][2][3],                  "RV_p5delchi2_lvl2"),
                (x->x[RVind][2][3][3],                  "RV_p5delchi2_lvl3"),

                (x->x[RVcom][1],                        "x_residuals_z_v0"),
                (x->x[RVcom][2],                        "x_residuals_v0"),
                # (x->x[RVcom][3],                        "x_skyLines_bright_v0"),
                (x->x[RVcom][3],                        "x_skyLines_faint_v0"),
                (x->x[RVcom][4],                        "x_skyContinuum_v0"),
                (x->x[RVcom][5],                        "x_starContinuum_v0"),
                # (x->x[RVcom][6],                        "x_starLineCor_v0"),
                (x->x[RVcom][6],                        "x_starLines_v0"),
                (x->x[RVcom][7],                        "x_starLineCof_v0"),
                (x->x[RVcom][8],                        "tot_p5chi2_v0"), 
                (x->x[RVcom][9],                        "apVisit_v0"),       
                                    
                (x->x[strpo],                           "x_starLines_err_v0"),    
            ]
            ## DIB Block
            DIBextract = []
            for (dibindx,dibw) in enumerate(dib_center_lst)
                dib = string(round(Int,dibw))
                dibind = string(dibindx)
                push!(DIBextract,[
                # Further chi2 refinement does not have fixed sizing because can hit grid edge
                (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][1][1]),        "DIB_pixoff_final_$(dibind)_$(dib)"),
                (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][1][2]),        "DIB_sigval_final_$(dibind)_$(dib)"),
                (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][3][1]),        "DIB_pixoff_disc_final_$(dibind)_$(dib)"),
                (x->Float64.(x[DIBind+dibsavesz*(dibindx-1)][1][3][2]),        "DIB_sigval_disc_final_$(dibind)_$(dib)"),
                (x->x[DIBind+dibsavesz*(dibindx-1)][1][2],                     "DIB_minchi2_final_$(dibind)_$(dib)"),
                (x->x[DIBind+dibsavesz*(dibindx-1)][1][6],                     "DIB_flag_$(dibind)_$(dib)"),
                (x->[x[DIBind+dibsavesz*(dibindx-1)][1][7:11]...],             "DIB_hess_var_$(dibind)_$(dib)"),
                                    
                (x->x[DIBind+dibsavesz*(dibindx-1)][2][1][3],                  "DIB_p5delchi2_lvl1_$(dibind)_$(dib)"),
                (x->x[DIBind+dibsavesz*(dibindx-1)][2][2][3],                  "DIB_p5delchi2_lvl2_$(dibind)_$(dib)"),
                (x->x[DIBind+dibsavesz*(dibindx-1)][2][3][3],                  "DIB_p5delchi2_lvl3_$(dibind)_$(dib)"),

                (x->x[EWind+dibsavesz*(dibindx-1)][1],                         "EW_dib_$(dibind)_$(dib)"),
                (x->x[EWind+dibsavesz*(dibindx-1)][2],                         "EW_dib_err_$(dibind)_$(dib)"),
                                    
                (x->x[DIBchi+dibsavesz*(dibindx-1)][1],                        "DIBchi2_residuals_$(dibind)_$(dib)"),

                (x->x[DIBcom+dibsavesz*(dibindx-1)][1],                        "x_residuals_z_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][2],                        "x_residuals_v1_$(dibind)_$(dib)"),
                # (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_bright_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_faint_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][4],                        "x_skyContinuum_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][5],                        "x_starContinuum_v1_$(dibind)_$(dib)"),
                # (x->x[DIBcom+dibsavesz*(dibindx-1)][6],                        "x_starLineCor_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][6],                        "x_starLines_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][7],                        "x_dib_v1_$(dibind)_$(dib)"),
                (x->x[DIBcom+dibsavesz*(dibindx-1)][8],                        "tot_p5chi2_v1_$(dibind)_$(dib)"),
                ])
            end
            extractlst = vcat(RVextract...,DIBextract...)
                
            hdr_dict = Dict(   
                    "pipeline"=>"apMADGICS.jl",
                    "git_branch"=>git_branch,   
                    "git_commit"=>git_commit,
            )           
            h5write(savename,"hdr","This is only a header")
            h5writeattr(savename,"hdr",hdr_dict)        
            for elelst in extractlst
                extractor(out,elelst[1],elelst[2],savename)
            end
        end
        GC.gc()
    end

    function extractor(x,elemap,elename,savename)
        len = length(x)
        exobj = elemap(x[1])
        outmat = zeros(eltype(exobj),size(exobj)...,len)
        for i=1:len
            flush(stdout)
            outmat[.. ,i] .= elemap(x[i])
        end
        h5write(savename,elename,outmat)
    end
end

batchsize = 10 #40
iterlst = []
Base.length(f::Iterators.Flatten) = sum(length, f.it)

for adjfibindx = 295:295 #1:600 #295, 245
    subiter = deserialize(prior_dir*"2024_01_19/outlists/dr17_dr17_star_input_lst_msked_"*lpad(adjfibindx,3,"0")*".jdat")
    subiterpart = Iterators.partition(subiter,batchsize)
    push!(iterlst,subiterpart)
end
ittot = Iterators.flatten(iterlst)
lenargs = length(ittot)
nwork = length(workers())
println("Batches to Do: $lenargs, number of workers: $nwork")
flush(stdout)

@showprogress pmap(multi_spectra_batch,ittot)
rmprocs(workers())

t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t0)); println("Total script runtime: $dt"); t_then = t_now; flush(stdout)