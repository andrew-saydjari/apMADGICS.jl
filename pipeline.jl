## This is the main pipeline that will batch over APOGEE files
# Author - Andrew Saydjari, CfA

import Pkg
Pkg.activate("./"); Pkg.instantiate(); Pkg.precompile()

using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager(launch_timeout=960.0))

activateout = @capture_out begin
    @everywhere begin
        import Pkg
        Pkg.activate("./")
    end
end

@everywhere begin
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, Interpolations, SparseArrays, ParallelDataTransfer, ThreadPinning, AstroTime
    prior_dir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch/working/"
    prior_dir2 = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/"
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
    
    using StatsBase, LinearAlgebra, ProgressMeter
    using BLISBLAS
    BLAS.set_num_threads(1)
end

# Task-Affinity CPU Locking in multinode SlurmContext
getinfo_worker(workerid::Int) = @getfrom workerid myid(), ThreadPinning.sched_getcpu(), gethostname()
idlst = getinfo_worker.(workers()); df = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]); push!(df,idlst...)
gdf = groupby(df,:hostname)
for sgdf in gdf, (sindx, sworker) in enumerate(sgdf.workerid)
    sendto(sworker, sindx=sindx)
    @spawnat sworker ThreadPinning.pinthread(sindx-1)
end
# Helpful Worker Info Printing
idlst = getinfo_worker.(workers()); df = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]); push!(df,idlst...)
gdf = groupby(df,:hostname); dfc = combine(gdf, nrow, :workerid => minimum, :workerid => maximum, :physcpu => minimum, :physcpu => maximum)
println("$(gethostname()) running Main on worker: $(myid()) cpu: $(ThreadPinning.sched_getcpu())")
for row in Tables.namedtupleiterator(dfc)
    println("$(row.hostname) running $(row.nrow) workers: $(row.workerid_minimum)->$(row.workerid_maximum) cpus: $(row.physcpu_minimum)->$(row.physcpu_maximum)")
end
flush(stdout)

println(BLAS.get_config())
flush(stdout)

using LibGit2
git_dir = src_dir
git_commit = LibGit2.head(git_dir)
git_repo = LibGit2.GitRepo(git_dir)
git_head = LibGit2.head(git_repo)
git_branch = LibGit2.shortname(git_head)
println("Running on branch: $git_branch, commit: $git_commit")
flush(stdout)

@passobj 1 workers() git_branch
@passobj 1 workers() git_commit

# These global allocations for the injest are messy... but we plan on changing the ingest
# relatively soon... so don't worry for now.
# SHOULD probably make a flag for using DD correction or not (do once we settle on the formalism)
@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg)
    
    c = 299792.458; # in km/s
    delLog = 6e-6; 
    pixscale = (10^(delLog)-1)*c;

    # nothing to do on size here, if anything expand
    f = h5open(prior_dir*"2023_03_07/precomp_dust_2_analyticDeriv.h5")
    global V_dib_noLSF = read(f["Vmat"])
    close(f)

    alpha = 1;
    f = h5open(prior_dir2*"2023_05_28/starLine_priors/APOGEE_stellar_kry_100_subpix_th22500.h5")
    global V_subpix_refLSF = alpha*read(f["Vmat"])
    close(f)

    # beta = 1;
    # f = h5open(prior_dir2*"2023_05_25/APOGEE_starCor_svd_50_subpix.h5")
    # global V_subpix_cor = beta*read(f["Vmat"])
    # global msk_starCor = convert.(Bool,read(f["msk_starCor"]))
    # close(f)

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

# it would be great to move this into a parameter file that is read for each run
@everywhere begin
    refine_iters = 1
    
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
    function pipeline_single_spectra(argtup; caching=true, cache_dir="../local_cache")
        ival = argtup[1]
        intup = argtup[2:end]
        out = []
        skycache = cache_skyname(intup,cache_dir=cache_dir)
        if (isfile(skycache) & caching)
            meanLocSky, VLocSky = deserialize(skycache)
        else
            meanLocSky, VLocSky = getSky4visit(intup)
            if caching
                dirName = splitdir(skycache)[1]
                if !ispath(dirName)
                    mkpath(dirName)
                end
                serialize(skycache,[meanLocSky, VLocSky])
            end
        end

        starcache = cache_starname(intup,cache_dir=cache_dir)
        if (isfile(starcache) & caching)
            fvec, fvarvec, cntvec, chipmidtimes = deserialize(starcache)
        else
            fvec, fvarvec, cntvec, chipmidtimes = stack_out(intup)
            if caching
                dirName = splitdir(starcache)[1]
                if !ispath(dirName)
                    mkpath(dirName)
                end
                serialize(starcache,[fvec, fvarvec, cntvec, chipmidtimes])
            end
        end
        framecnts = maximum(cntvec)
        simplemsk = (cntvec.==framecnts) .& skymsk;
        
        starscale = if count(simplemsk .& (.!isnan.(fvec)))==0
            NaN
        else
            abs(nanmedian(fvec[simplemsk]))
        end
        push!(out,(count(simplemsk), starscale, framecnts, chipmidtimes)) # 1

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
        V_skyline_c = V_skyline
        V_skyline_r = V_skyline_c[simplemsk,:]
        V_locSky_c = VLocSky
        V_locSky_r = V_locSky_c[simplemsk,:]
        V_starCont_c = starscale*V_starcont
        V_starCont_r = V_starCont_c[simplemsk,:]

        ## Solve RV of Star
        # compute stellar continuum to modify stellar line prior
        Vcomb_cur = hcat(V_skyline_r,V_locSky_r,V_starCont_r);
        Ctotinv_cur = LowRankMultMatIP([Ainv,Vcomb_cur],wood_precomp_mult_mat([Ainv,Vcomb_cur],(size(Ainv,1),size(V_subpix_comb,2))),wood_fxn_mult,wood_fxn_mult_mat!);
        x_comp_lst = deblend_components_all(Ctotinv_cur, Xd_obs, (V_starCont_r,))

        starCont_Mscale = x_comp_lst[1]
        pre_Vslice = zeros(count(simplemsk),size(V_subpix_comb,2))
        chi2_wrapper_partial = Base.Fix2(chi2_wrapper,(simplemsk,Ctotinv_cur,Xd_obs,starCont_Mscale,V_subpix_comb,pre_Vslice))
        lout = sampler_1d_hierarchy_var(chi2_wrapper_partial,slvl_tuple,minres=1//10,stepx=1)
        push!(out,lout) # 2

        # update the Ctotinv to include the stellar line component (iterate to refine starCont_Mscale)
        svalc = lout[1][3]
        for i=1:refine_iters
            Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_cur.matList[1],simplemsk,starCont_Mscale,Vcomb_cur,V_subpix,V_subpix_refLSF)
            x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,))
            starCont_Mscale = x_comp_lst[1]
        end
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines_asym(svalc,Ctotinv_cur.matList[1],simplemsk,starCont_Mscale,Vcomb_cur,V_subpix,V_subpix_refLSF)
        
        # do a component save without the 15273 DIB
        x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
            (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_r),
            (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_c),
        )
        push!(out,x_comp_lst[1]'*(Ainv*x_comp_lst[1])) # 3
        x_comp_out = [nanify(x_comp_lst[1],simplemsk)./sqrt.(fvarvec), nanify(x_comp_lst[1],simplemsk), nanify(x_comp_lst[2],simplemsk), 
                        nanify(x_comp_lst[3].+meanLocSky[simplemsk],simplemsk), nanify(x_comp_lst[4],simplemsk),
                        x_comp_lst[5:end]...]
        push!(out,x_comp_out) # 4
        dflux_starlines = sqrt_nan.(get_diag_posterior_from_prior_asym(Ctotinv_fut, V_starlines_c, V_starlines_r))
        push!(out,dflux_starlines) # 5
                
        # prepare multiplicative factors for DIB prior
        x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,V_starlines_r))
        starCont_Mscale = x_comp_lst[1]
        starFull_Mscale = x_comp_lst[1].+x_comp_lst[2]
        
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines2_asym(svalc,Ctotinv_cur.matList[1],simplemsk,starCont_Mscale,Vcomb_cur,V_subpix,V_subpix_refLSF)
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
            lout = sampler_2d_hierarchy_var(chi2_wrapper_partial,lvltuple_dib)
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
                opt_tup,Ctotinv_cur.matList[1],simplemsk,starFull_Mscale,Vcomb_cur,V_dib,V_dib_noLSF,scan_offset)

            x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
                (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_r, V_dibr),
                (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_c, V_dibc),
            )
            push!(out,x_comp_lst[1]'*(Ainv*x_comp_lst[1])) # 8
            # I am not sure that during production we really want to run and output full sets of components per DIB
            # I would like to fill NaNs in chip gaps for the sky/continuum components
            # revisit that when we revisit the interpolations before making other fiber priors
            x_comp_out = [nanify(x_comp_lst[1],simplemsk)./sqrt.(fvarvec), nanify(x_comp_lst[1],simplemsk), nanify(x_comp_lst[2],simplemsk), 
                        nanify(x_comp_lst[3].+meanLocSky[simplemsk],simplemsk), nanify(x_comp_lst[4],simplemsk),
                        x_comp_lst[5:end]...]

            push!(out,x_comp_out) # 9
        end
                        
        return out
    end
end

@everywhere begin
    function multi_spectra_batch(indsubset; out_dir="../outdir")
        ### Set up
        out = []
        startind = indsubset[1][1]
        tele = indsubset[1][2]
        fiberindx = indsubset[1][end]
        teleind = (tele == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiberindx

        ### Save and cache restart handling
        savename = join([out_dir,lpad(adjfibindx,3,"0"),"apMADGICS_fiber_"*lpad(adjfibindx,3,"0")*"_batch_"*lpad(startind,7,"0")*".h5"],"/")
        dirName = splitdir(savename)[1]
        if !ispath(dirName)
            mkpath(dirName)
        end
        if !isfile(savename)
            ### Need to load the priors here
            f = h5open(prior_dir2*"2023_06_01/sky_priors/APOGEE_skycont_svd_30_f"*lpad(adjfibindx,3,"0")*".h5")
            global V_skycont = read(f["Vmat"])
            chebmsk_exp = convert.(Bool,read(f["chebmsk_exp"]))
            close(f)

            f = h5open(prior_dir2*"2023_06_01/sky_priors/APOGEE_skyline_svd_120_f"*lpad(adjfibindx,3,"0")*".h5")
            global V_skyline = read(f["Vmat"])
            submsk = convert.(Bool,read(f["submsk"]))
            close(f)

            global skymsk = chebmsk_exp .& submsk #.& msk_starCor;

            f = h5open(prior_dir2*"2023_06_01/star_priors/APOGEE_starcont_svd_60_f"*lpad(adjfibindx,3,"0")*".h5")
            global V_starcont = read(f["Vmat"])
            close(f)

            # can consider changing dimension at the full DR17 reduction stage
            f = h5open(prior_dir2*"2023_05_28/starLine_priors/APOGEE_stellar_kry_100_subpix_"*lpad(adjfibindx,3,"0")*".h5")
            global V_subpix = alpha*read(f["Vmat"])
            close(f)
            # global V_subpix_comb = hcat(V_subpix,V_subpix_cor)
            global V_subpix_comb = V_subpix

            f = h5open(prior_dir*"2023_04_03/dib_priors/precomp_dust_2_analyticDerivLSF_"*lpad(adjfibindx,3,"0")*".h5")
            global V_dib = read(f["Vmat"])
            close(f)

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
                (x->adjfibindx,                         "adjfiberindx"),

                (x->Float64.(x[RVind][1][1]),           "RV_pixoff_final"),
                (x->x[RVind][1][2],                     "RV_minchi2_final"),
                (x->x[RVind][1][6],                     "RV_flag"),
                (x->x[RVind][1][7],                     "RV_pix_var"),
                                    
                (x->x[RVchi][1],                        "RVchi2_residuals"),
                                    
                (x->x[RVind][2][1][3],                  "RV_p5delchi2_lvl1"),
                (x->x[RVind][2][2][3],                  "RV_p5delchi2_lvl2"),
                (x->x[RVind][2][3][3],                  "RV_p5delchi2_lvl3"),

                (x->x[RVcom][1],                        "x_residuals_z_v0"),
                (x->x[RVcom][2],                        "x_residuals_v0"),
                (x->x[RVcom][3],                        "x_skyLines_v0"),
                (x->x[RVcom][4],                        "x_skyContinuum_v0"),
                (x->x[RVcom][5],                        "x_starContinuum_v0"),
                # (x->x[RVcom][6],                        "x_starLineCor_v0"),
                (x->x[RVcom][6],                        "x_starLines_v0"),
                (x->x[RVcom][7],                        "tot_p5chi2_v0"),       
                                    
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
                (x->x[DIBcom+dibsavesz*(dibindx-1)][3],                        "x_skyLines_v1_$(dibind)_$(dib)"),
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
for adjfibindx=295:295
    subiter = deserialize(prior_dir*"2023_04_04/star_input_lists/star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat")
    subiterpart = Iterators.partition(subiter,batchsize)
    push!(iterlst,subiterpart)
end
ittot = Iterators.flatten(iterlst)
lenargs = length(ittot)
nwork = length(workers())
println("Batches to Do: $lenargs, number of workers: $nwork")
flush(stdout)

@showprogress pmap(multi_spectra_batch,ittot)
