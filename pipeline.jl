## This is the main pipeline that will batch over APOGEE files
import Pkg
Pkg.activate("./"); Pkg.instantiate(); Pkg.precompile()

using Distributed, SlurmClusterManager, LibGit2
addprocs(SlurmManager())
        
@everywhere println("hello from $(myid()):$(gethostname())")
flush(stdout)

@everywhere begin
    import Pkg
    Pkg.activate("./")
end
flush(stdout)

@everywhere begin
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays, Interpolations, SparseArrays
    include("src/utils.jl")
    include("src/gridSearch.jl")
    include("src/componentAndPosteriors.jl")
    include("src/fileNameHandling.jl")
    include("src/ingest.jl")
    include("src/lowRankPrescription.jl")
    include("src/marginalizeEW.jl")
    include("src/spectraInterpolation.jl")
    include("src/chi2Wrappers.jl")
    
    using StatsBase, LinearAlgebra, ProgressMeter
    BLAS.set_num_threads(1)
end

git_dir = "./"
git_commit = LibGit2.head(git_dir)
git_repo = LibGit2.GitRepo(git_dir)
git_head = LibGit2.head(git_repo)
git_branch = LibGit2.shortname(git_head)
println("Running on branch: $git_branch, commit: $git_commit")
flush(stdout)

# These global allocations for the injest are messy... but we plan on changing the ingest
# relatively soon... so don't worry for now.
@everywhere begin
    wavetarg = 10 .^range(start=(4.179-125*6.0e-6),step=6.0e-6,length=8575+125);
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
    f = h5open("../../2023_02_28/APOGEE_skycont_svd_150_f295.h5")
    V_skycont = f["Vmat"][:,1:30]
    close(f)

    # pretty happy here, could be convinced to decrease a little bit
    f = h5open("../../2023_02_28/APOGEE_skyline_kry_150_f295.h5")
    V_skyline = f["Vmat"][:,1:100]
    close(f)

    f = h5open("../../2023_03_03/APOGEE_starcont_svd_150_f295.h5")
    V_starcont = f["Vmat"][:,1:60]
    close(f)

    # hard to test and decide to decrease without doing a batch over a large range of stellar types
    # can consider dropping at the full fiber reduction stage
    f = h5open("../../2023_03_06/APOGEE_stellar_svd_50_f295_lite_subpix_zerocent.h5")
    V_subpix = read(f["Vmat"])
    close(f)

    # nothing to do on size here, if anything expand
    f = h5open("../../2023_03_07/precomp_dust_2_analyticDeriv.h5")
    V_dib = read(f["Vmat"])
    close(f)
end

# it would be great to move this into a parameter file that is red for each run
@everywhere begin
    # Star Wave
    lvl1 = -70:1//2:70
    lvl2 = -8:2//10:8
    lvl3 = -3:1//10:3
    slvl_tuple = (lvl1,lvl2,lvl3)
    # tuple1dprint(slvl_tuple)

    # (Wave, Sig) DIB
    dib_center = 15273
    lvl1d = ((-150:4:150),(18//10:18//10))
    lvl2d = ((0:0), (-7//5:4//100:11/5))
    lvl3d = ((-18:2//10:18), (0:0))
    lvl4d = ((0:0), (-90//100:2//100:90//100))
    lvl5d = ((-4//10:1//10:4//10), (-4//100:1//100:4//100));
    lvltuple = (lvl1d, lvl2d, lvl3d, lvl4d, lvl5d);
    # tuple2dprint(lvltuple)

    # Flux marginalize region
    sigMarg0 = -50//100:10//100:50//100
    svalMarg0 = -0//10:1//10:0//10;
end

@everywhere begin
    function pipeline_single_spectra(argtup; caching=true)
        ival = argtup[1]
        intup = argtup[2:end]
        out = []
        skycache = cache_skyname(intup)
        if (isfile(skycache) & caching)
            meanLocSky, VLocSky = deserialize(skycache)
        else
            try
                meanLocSky, VLocSky = getSky4visit(intup)
                if caching
                    serialize(skycache,[meanLocSky, VLocSky])
                end
            catch
                println(intup)
            end
        end

        starcache = cache_starname(intup)
        if (isfile(starcache) & caching)
            fvec, fvarvec, cntvec = deserialize(starcache)
        else
            fvec, fvarvec, cntvec = stack_out(intup)
            if caching
                serialize(starcache,[fvec, fvarvec, cntvec])
            end
        end
        simplemsk = (cntvec.==maximum(cntvec));
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
        Vcomb = hcat(V_skyline_r,V_locSky_r,V_starCont_r);
        Ctotinv = LowRankMultMatIP([Ainv,Vcomb],wood_precomp_mult_mat([Ainv,Vcomb],(size(Ainv,1),size(V_subpix,2))),wood_fxn_mult,wood_fxn_mult_mat!);
        x_comp_lst = deblend_components_all(Ctotinv, Xd_obs, (V_starCont_r,))

        starCont_Mscale = Diagonal(x_comp_lst[1])
        chi2_wrapper_partial = Base.Fix2(chi2_wrapper,(simplemsk,Ctotinv,Xd_obs,starCont_Mscale,V_subpix))
        lout = sampler_1d_hierarchy_var(chi2_wrapper_partial,slvl_tuple,minres=1//10,stepx=1)
        push!(out,lout)

        # update the Ctot inv to include the stellar line component (could iterate to refine starContScale)
        svalc = lout[1][3]
        Ctotinv, Vcomb, V_starlines_c, V_starlines_r = update_Ctotinv_Vstarstarlines(svalc,Ctotinv.matList[1],simplemsk,starCont_Mscale,Vcomb,V_subpix)
        x_comp_lst = deblend_components_all(Ctotinv, Xd_obs, (V_starCont_r,V_starlines_r))
        starCont_Mscale = Diagonal(x_comp_lst[1])
        starFull_Mscale = Diagonal(x_comp_lst[1].+x_comp_lst[2])

        ## Solve DIB parameters (for just 15273)
        # one of the main questions is how many time to compute components and where
        chi2_wrapper_partial = Base.Fix2(chi2_wrapper2d,(simplemsk,Ctotinv,Xd_obs,wave_obs,starFull_Mscale,Vcomb,V_dib,dib_center))
        lout = sampler_2d_hierarchy_var(chi2_wrapper_partial,lvltuple)
        opt_tup = lout[1][3]
        push!(out,lout)

        ## Shift the marginalization sampling (should this be wrapped inside the function?)
        # especially because we need to do bounds handling
        svalMarg = svalMarg0 .+ lout[1][3][1]
        sigMarg = shift_trim_range(sigMarg0,lout[1][3][2]; minv=4//10, maxv=4)
        samp_lst = Iterators.product(svalMarg,sigMarg)

        intupf = (simplemsk,Ctotinv,Xd_obs,wave_obs,starFull_Mscale,Vcomb,V_dib,dib_center)
        chi2lst, fluxlst, dfluxlst = sample_chi2_flux_dflux(samp_lst,intupf)
        refchi2val = minimum(chi2lst) #this should just be set to the min found at the 2d step
        lout = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
        push!(out,lout)

        # Compute some final components for export
        Ctotinv, Vcomb, V_dibc, V_dibr = update_Ctotinv_Vdib(
            opt_tup,Ctotinv.matList[1],simplemsk,starFull_Mscale,Vcomb,V_dib)

        x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, 
            (A, V_skyline_r, V_locSky_r, V_starCont_r, V_starlines_r, V_dibr),
            (A, V_skyline_c, V_locSky_c, V_starCont_c, V_starlines_c, V_dibc),
        )

        # I would like to fill NaNs in chip gaps for the sky/continuum components
        # revisit that when we revisit the interpolations before making other fiber priors
        x_comp_out = [nanify(x_comp_lst[1],simplemsk), x_comp_lst[2], x_comp_lst[3].+meanLocSky, x_comp_lst[4:end]...]

        push!(out,x_comp_out)
#         push!(out,(wave_obs,fvarvec[simplemsk],simplemsk))
#         push!(out,(meanLocSky, VLocSky))
        return out
    end

    function multi_spectra_batch(indsubset; fibnum=295, out_dir="../outdir/")
        out = []
        for (ind,indval) in enumerate(indsubset)
            push!(out,pipeline_single_spectra(indval; caching=true))
        end

        startind = indsubset[1][1]
        savename = out_dir*"apMADGICS_fiber_"*lpad(fibnum,3,"0")*"_batch_"*lpad(startind,7,"0")*".h5"

        extractlst = [
            (x->x[1][2][1][3],                  "RV_p5delchi2_lvl1"),
            (x->x[1][2][2][3],                  "RV_p5delchi2_lvl2"),
            (x->x[1][2][3][3],                  "RV_p5delchi2_lvl3"),

            (x->x[1][1][1],                     "RV_minchi2_final"),
            (x->x[1][1][2],                     "RV_pixoff_final"),
            (x->x[1][1][5],                     "RV_flag"),
            (x->x[1][1][6],                     "RV_pix_var"),

            (x->x[2][2][1][3],                  "DIB_p5delchi2_lvl1"),
            (x->x[2][2][2][3],                  "DIB_p5delchi2_lvl2"),
            (x->x[2][2][3][3],                  "DIB_p5delchi2_lvl3"),
            # These do not have fixed sizing because they can hit the grid edge for sigma... need to ponder if/how to handle
#             (x->x[2][2][4][3],      "DIB_p5delchi2_lvl4"),
#             (x->x[2][2][5][3],      "DIB_p5delchi2_lvl5"),

            (x->x[2][1][1],                     "DIB_minchi2_final"),
            (x->Float64.(x[2][1][2][1]),        "DIB_pixoff_final"),
            (x->Float64.(x[2][1][2][2]),        "DIB_sigval_final"),
            (x->x[2][1][5],                     "DIB_flag"),
            (x->[x[2][1][6:10]...],             "DIB_hess_var"),

            (x->x[3][1],                        "EW_dib"),
            (x->x[3][2],                        "EW_dib_err"),

            (x->x[4][1],                        "x_residuals_v1"),
            (x->x[4][2],                        "x_skyLines_v1"),
            (x->x[4][3],                        "x_skyContinuum_v1"),
            (x->x[4][4],                        "x_starContinuum_v1"),
            (x->x[4][5],                        "x_starLines_v1"),
            (x->x[4][6],                        "x_dib_v1"),
        ]

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
