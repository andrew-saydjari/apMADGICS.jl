## This is the main pipeline that will batch over APOGEE files
# Author - Andrew Saydjari, CfA

import Pkg
Pkg.activate("./"); Pkg.instantiate(); Pkg.precompile()

using Distributed, SlurmClusterManager
addprocs(SlurmManager(launch_timeout=960.0))
        
@everywhere println("hello from $(myid()):$(gethostname())")
flush(stdout)

@everywhere begin
    import Pkg
    Pkg.activate("./")
end
flush(stdout)

@everywhere begin
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, Interpolations, SparseArrays, ParallelDataTransfer
    prior_dir = "../../"
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
@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg)
    
    c = 299792.458; # in km/s
    delLog = 6e-6; 
    pixscale = (10^(delLog)-1)*c;

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

# This overhead is going to depend on fiber number soon, so this will move inside the multispectra wrapper
@everywhere begin
    # pretty happy at here, revisit if we incoporate tellurics more consistently
    f = h5open(prior_dir*"2023_03_29/APOGEE_skycont_svd_150_f295.h5")
    V_skycont = f["Vmat"][:,1:30]
    chebmsk_exp = convert.(Bool,read(f["chebmsk_exp"]))
    close(f)

    # pretty happy here, could be convinced to decrease a little bit
    f = h5open(prior_dir*"2023_03_29/APOGEE_skyline_svd_150_f295.h5")
    V_skyline = f["Vmat"][:,1:100]
    submsk = convert.(Bool,read(f["submsk"]))
    close(f)
    
    skymsk = chebmsk_exp .& submsk;

    f = h5open(prior_dir*"2023_03_29/APOGEE_starcont_svd_150_f295.h5")
    V_starcont = f["Vmat"][:,1:60]
    close(f)

    # hard to test and decide to decrease without doing a batch over a large range of stellar types
    # can consider dropping at the full fiber reduction stage
    f = h5open(prior_dir*"2023_03_06/APOGEE_stellar_svd_50_f295_lite_subpix_zerocent.h5")
    V_subpix = read(f["Vmat"])
    close(f)

    # nothing to do on size here, if anything expand
    f = h5open(prior_dir*"2023_03_07/precomp_dust_2_analyticDeriv.h5")
    V_dib_noLSF = read(f["Vmat"])
    close(f)
        
    f = h5open(prior_dir*"2023_03_23/precomp_dust_2_analyticDerivLSF.h5")
    V_dib = read(f["Vmat"])
    close(f)
end

# it would be great to move this into a parameter file that is red for each run
@everywhere begin
    refine_iters = 1
    
    # Star Wave
    lvl1 = -70:1//2:70
    lvl2 = -8:2//10:8
    lvl3 = -3:1//10:3
    slvl_tuple = (lvl1,lvl2,lvl3)
    # tuple1dprint(slvl_tuple)

    # (Wave, Sig) DIB
    dib_center_lst = [15273]#, 15653]
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
#             try
                meanLocSky, VLocSky = getSky4visit(intup)
#                 if caching
#                     serialize(skycache,[meanLocSky, VLocSky])
#                 end
#             catch
#                 println(intup)
#             end
        end

        starcache = cache_starname(intup,cache_dir=cache_dir)
        if (isfile(starcache) & caching)
            fvec, fvarvec, cntvec = deserialize(starcache)
        else
            fvec, fvarvec, cntvec = stack_out(intup)
            if caching
                serialize(starcache,[fvec, fvarvec, cntvec])
            end
        end
        simplemsk = (cntvec.==maximum(cntvec)) .& skymsk;
        fvec./=maximum(cntvec)
        fvarvec./=(maximum(cntvec)^2)
        
        starscale = abs(nanmedian(fvec[simplemsk]))

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
        Ctotinv_cur = LowRankMultMatIP([Ainv,Vcomb_cur],wood_precomp_mult_mat([Ainv,Vcomb_cur],(size(Ainv,1),size(V_subpix,2))),wood_fxn_mult,wood_fxn_mult_mat!);
        x_comp_lst = deblend_components_all(Ctotinv_cur, Xd_obs, (V_starCont_r,))

        starCont_Mscale = x_comp_lst[1]
        pre_Vslice = zeros(count(simplemsk),size(V_subpix,2))
        chi2_wrapper_partial = Base.Fix2(chi2_wrapper,(simplemsk,Ctotinv_cur,Xd_obs,starCont_Mscale,V_subpix,pre_Vslice))
        lout = sampler_1d_hierarchy_var(chi2_wrapper_partial,slvl_tuple,minres=1//10,stepx=1)
        push!(out,lout) # 1

        # update the Ctotinv to include the stellar line component (iterate to refine starCont_Mscale)
        svalc = lout[1][3]
        for i=1:refine_iters
            Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines(svalc,Ctotinv_cur.matList[1],simplemsk,starCont_Mscale,Vcomb_cur,V_subpix)
            x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,))
            starCont_Mscale = x_comp_lst[1]
        end
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines(svalc,Ctotinv_cur.matList[1],simplemsk,starCont_Mscale,Vcomb_cur,V_subpix)
        
        # do a component save without the 15273 DIB
        x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
            (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_r),
            (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_c),
        )
        push!(out,x_comp_lst[1]'*(Ainv*x_comp_lst[1])) # 2
        x_comp_out = [nanify(x_comp_lst[1],simplemsk), nanify(x_comp_lst[2],simplemsk), 
                nanify(x_comp_lst[3].+meanLocSky[simplemsk],simplemsk), nanify(x_comp_lst[4],simplemsk),
                x_comp_lst[5]]
        push!(out,x_comp_out) # 3
        dflux_starlines = sqrt_nan.(get_diag_posterior_from_prior_asym(Ctotinv_fut, V_starlines_c, V_starlines_r))
        push!(out,dflux_starlines) # 4
                
        # prepare multiplicative factors for DIB prior
        x_comp_lst = deblend_components_all(Ctotinv_fut, Xd_obs, (V_starCont_r,V_starlines_r))
        starCont_Mscale = x_comp_lst[1]
        starFull_Mscale = x_comp_lst[1].+x_comp_lst[2]
        
        Ctotinv_fut, Vcomb_fut, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines(svalc,Ctotinv_cur.matList[1],simplemsk,starCont_Mscale,Vcomb_cur,V_subpix)
        Ctotinv_cur, Ctotinv_fut = Ctotinv_fut, Ctotinv_cur; Vcomb_cur, Vcomb_fut = Vcomb_fut, Vcomb_cur # swap to updated covariance finally
        
        # currently, this is modeling each DIB seperately... I think we want to change this later, just easier parallel structure
        pre_Vslice = zeros(count(simplemsk),size(V_dib,2))
        for dib_ind = 1:length(dib_center_lst) # eventually need to decide if these are cumulative or not
            dib_center = dib_center_lst[dib_ind]
            scan_offset = findmin(abs.(wavetarg.-dib_center_lst[dib_ind]))[2].-findmin(abs.(wavetarg.-dib_center_lst[1]))[2]
            
            ## Solve DIB parameters (for just 15273), not any more, just a single DIB
            # one of the main questions is how many time to compute components and where
            chi2_wrapper_partial = Base.Fix2(chi2_wrapper2d,(simplemsk,Ctotinv_cur,Xd_obs,wave_obs,starFull_Mscale,Vcomb_cur,V_dib,pre_Vslice,dib_center,scan_offset))
            lout = sampler_2d_hierarchy_var(chi2_wrapper_partial,lvltuple)
            opt_tup = lout[1][3]
            push!(out,lout) # 5, 9

            ## Shift the marginalization sampling (should this be wrapped inside the function?)
            # especially because we need to do bounds handling
            svalMarg = svalMarg0 .+ opt_tup[1]
            sigMarg = shift_trim_range(sigMarg0,opt_tup[2]; minv=4//10, maxv=4)
            samp_lst = Iterators.product(svalMarg,sigMarg)

            intupf = (simplemsk,Ctotinv_cur,Xd_obs,wave_obs,starFull_Mscale,Vcomb_cur,V_dib,pre_Vslice,dib_center,scan_offset)
            chi2lst, fluxlst, dfluxlst = sample_chi2_flux_dflux(samp_lst,intupf) #shouldn't this take chi2_wrapper_partial as an argument?
            refchi2val = minimum(chi2lst) #this should just be set to the min found at the 2d step
            lout = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
            push!(out,lout) # 6, 10

            # Compute some final components for export (still need to implement DIB iterative refinement)
            Ctotinv_fut, Vcomb_fut, V_dibc, V_dibr = update_Ctotinv_Vdib_asym(
                opt_tup,Ctotinv_cur.matList[1],simplemsk,starFull_Mscale,Vcomb_cur,V_dib,V_dib_noLSF,scan_offset)

            x_comp_lst = deblend_components_all_asym_tot(Ctotinv_fut, Xd_obs, 
                (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_r, V_dibr),
                (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_c, V_dibc),
            )
            push!(out,x_comp_lst[1]'*(Ainv*x_comp_lst[1])) # 7, 11
            # I am not sure that during production we really want to run and output full sets of components per DIB
            # I would like to fill NaNs in chip gaps for the sky/continuum components
            # revisit that when we revisit the interpolations before making other fiber priors
            x_comp_out = [nanify(x_comp_lst[1],simplemsk), nanify(x_comp_lst[2],simplemsk), 
                        nanify(x_comp_lst[3].+meanLocSky[simplemsk],simplemsk), nanify(x_comp_lst[4],simplemsk),
                        x_comp_lst[5:end]...]

            push!(out,x_comp_out) # 8, 12
        end
                        
        push!(out,count(simplemsk)) # 13
#         push!(out,(wave_obs,fvarvec[simplemsk],simplemsk))
#         push!(out,(meanLocSky, VLocSky))
        return out
    end
end
@everywhere begin
    function multi_spectra_batch(indsubset; fibnum=295, out_dir="../outdir/")
        out = []
        for (ind,indval) in enumerate(indsubset)
            push!(out,pipeline_single_spectra(indval; caching=true))
        end

        startind = indsubset[1][1]
        savename = out_dir*"apMADGICS_fiber_"*lpad(fibnum,3,"0")*"_batch_"*lpad(startind,7,"0")*".h5"
        
        RVind = 1
        RVchi = 2
        RVcom = 3
        strpo = 4
            
        DIBin1 = 5
        EWind1 = 6
        DIBch1 = 7
        DIBco1 = 8
            
#         DIBin2 = 9
#         EWind2 = 10
#         DIBch2 = 11
#         DIBco2 = 12
            
        metai = 9 #13
        extractlst = [
            (x->x[RVind][1][1],                     "RV_pixoff_final"),
            (x->x[RVind][1][2],                     "RV_minchi2_final"),
            (x->x[RVind][1][6],                     "RV_flag"),
            (x->x[RVind][1][7],                     "RV_pix_var"),
                                
            (x->x[RVchi][1],                        "RVchi2_residuals"),
                                
            (x->x[RVind][2][1][3],                  "RV_p5delchi2_lvl1"),
            (x->x[RVind][2][2][3],                  "RV_p5delchi2_lvl2"),
            (x->x[RVind][2][3][3],                  "RV_p5delchi2_lvl3"),

            (x->x[RVcom][1],                        "x_residuals_v0"),
            (x->x[RVcom][2],                        "x_skyLines_v0"),
            (x->x[RVcom][3],                        "x_skyContinuum_v0"),
            (x->x[RVcom][4],                        "x_starContinuum_v0"),
            (x->x[RVcom][5],                        "x_starLines_v0"),
            (x->x[RVcom][6],                        "tot_p5chi2_v0"),       
                                
            (x->x[strpo],                           "x_starLines_err_v0"),    
            
            # we should automate this populating the name from DIB center wave
            (x->Float64.(x[DIBin1][1][1][1]),        "DIB_pixoff_final_15273"),
            (x->Float64.(x[DIBin1][1][1][2]),        "DIB_sigval_final_15273"),
            (x->x[DIBin1][1][2],                     "DIB_minchi2_final_15273"),
            (x->x[DIBin1][1][6],                     "DIB_flag_15273"),
            (x->[x[DIBin1][1][7:11]...],             "DIB_hess_var_15273"),
                                
            (x->x[DIBin1][2][1][3],                  "DIB_p5delchi2_lvl1_15273"),
            (x->x[DIBin1][2][2][3],                  "DIB_p5delchi2_lvl2_15273"),
            (x->x[DIBin1][2][3][3],                  "DIB_p5delchi2_lvl3_15273"),
            # These do not have fixed sizing because they can hit the grid edge for sigma... need to ponder if/how to handle
#             (x->x[DIBin][2][4][3],      "DIB_p5delchi2_lvl4"),
#             (x->x[DIBin][2][5][3],      "DIB_p5delchi2_lvl5"),

            (x->x[EWind1][1],                        "EW_dib_15273"),
            (x->x[EWind1][2],                        "EW_dib_err_15273"),
                                
            (x->x[DIBch1][1],                        "DIBchi2_residuals_15273"),

            (x->x[DIBco1][1],                        "x_residuals_v1_15273"),
            (x->x[DIBco1][2],                        "x_skyLines_v1_15273"),
            (x->x[DIBco1][3],                        "x_skyContinuum_v1_15273"),
            (x->x[DIBco1][4],                        "x_starContinuum_v1_15273"),
            (x->x[DIBco1][5],                        "x_starLines_v1_15273"),
            (x->x[DIBco1][6],                        "x_dib_v1_15273"),
            (x->x[DIBco1][7],                        "tot_p5chi2_v1_15273"),
                
#             (x->Float64.(x[DIBin2][1][1][1]),        "DIB_pixoff_final_15673"),
#             (x->Float64.(x[DIBin2][1][1][2]),        "DIB_sigval_final_15673"),
#             (x->x[DIBin2][1][2],                     "DIB_minchi2_final_15673"),
#             (x->x[DIBin2][1][6],                     "DIB_flag_15673"),
#             (x->[x[DIBin2][1][7:11]...],             "DIB_hess_var_15673"),
                                
#             (x->x[DIBin2][2][1][3],                  "DIB_p5delchi2_lvl1_15673"),
#             (x->x[DIBin2][2][2][3],                  "DIB_p5delchi2_lvl2_15673"),
#             (x->x[DIBin2][2][3][3],                  "DIB_p5delchi2_lvl3_15673"),

#             (x->x[EWind2][1],                        "EW_dib_15673"),
#             (x->x[EWind2][2],                        "EW_dib_err_15673"),
                                
#             (x->x[DIBch2][1],                        "DIBchi2_residuals_15673"),

#             (x->x[DIBco2][1],                        "x_residuals_v1_15673"),
#             (x->x[DIBco2][2],                        "x_skyLines_v1_15673"),
#             (x->x[DIBco2][3],                        "x_skyContinuum_v1_15673"),
#             (x->x[DIBco2][4],                        "x_starContinuum_v1_15673"),
#             (x->x[DIBco2][5],                        "x_starLines_v1_15673"),
#             (x->x[DIBco2][6],                        "x_dib_v1_15673"),
#             (x->x[DIBco2][7],                        "tot_p5chi2_v1_15673"),
                                
            (x->x[metai],                           "data_pix_cnt"), #this is a DOF proxy, but I think our more careful info/pixel analysis would be better
        ]
             
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

input_list = deserialize("../input_list.jl")
itarg = Iterators.partition(input_list,10)
larg = length(itarg)
nwork = length(workers())
println("Batches to Do: $larg, number of workers: $nwork")
flush(stdout)

@showprogress pmap(multi_spectra_batch,itarg)
