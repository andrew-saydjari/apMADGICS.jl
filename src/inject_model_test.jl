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
    correlated_noise = true

    dib_center_lambda_lst = [15273] #,15672]
    dib_ew_range = (-1.5,0)
    dib_sig_range = (0.7,3.7)
    dib_vel_range = (-450, 450) # km/s

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    prior_dict["out_dir"] = prior_dir*"2024_03_11/inject_15273only_295_corr/"
    prior_dict["inject_cache_dir"] = prior_dir*"2024_03_11/inject_local_cache_15273only_295_corr/"
    prior_dict["local_cache"] = prior_dir*"2024_03_11/local_cache_inject_corr/"

    prior_dict["past_run"] = prior_dir*"2024_02_21/outdir_wu_295_th/apMADGICS_out.h5" # used for StarScale distribution only
    prior_dict["korg_run_path"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/starLine_disk_KnackedKorg/"

    prior_dict["sqrt_corr_kernels"] = prior_dir*"2024_03_11/sqrt_corr_kernel_"

    # Sky Sources
    prior_dict["sky_runlist"] = prior_dir*"2024_02_21/outlists/sky/dr17_dr17_sky_input_lst_plate_msked_"
    prior_dict["starContSamples"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/tell_prior_disk/starCont_"
    prior_dict["skyContSamples"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/sky_prior_disk/skycont_"
    prior_dict["chebmsk_exp"] = prior_dir*"2024_02_21/apMADGICS.jl/src/prior_build/sky_prior_disk/chebmsk_exp_"

    prior_dict["LSF_mat_APO"] = prior_dir0*"2023_04_01/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_01 by AKS
    prior_dict["LSF_mat_LCO"] = prior_dir0*"2023_04_07/mat_lsf_out/sp_combolsfmat_norm_" # last made 2023_04_07 by AKS

    prior_dict["chip_fluxdep_err_correction"] = src_dir*"data/chip_fluxdep_err_correction.jdat"

    rnd_seed = 695
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
    savename = prior_dict["starContSamples"]*lpad(adjfibindx,3,"0")*".jdat"
    starContSamples = deserialize(savename);

    savename = prior_dict["chebmsk_exp"]*lpad(adjfibindx,3,"0")*".jdat"
    chebmsk_exp = deserialize(savename);

    adjfiberindx = reader(prior_dict["past_run"],"adjfiberindx");
    mskFIB = (adjfiberindx.==295)
    starscale = reader(prior_dict["past_run"],"starscale")
    starscale_red = starscale[..,mskFIB];

    savename = prior_dict["skyContSamples"]*lpad(adjfibindx,3,"0")*".jdat"
    skyContSamples_raw = deserialize(savename)
    msk_skyCont = map(x->(!any(isnan.(x)))&(!any(x.<0)),eachcol(skyContSamples_raw));
    skyContSamples = skyContSamples_raw[:,msk_skyCont]

    err_correct_Dict = deserialize(prior_dict["chip_fluxdep_err_correction"])

    sqrtMats = Dict{String,Matrix}()
    sqrtMats["apo"] = deserialize(prior_dict["sqrt_corr_kernels"]*"APO.jdat")
    # sqrtMats["lco"] = deserialize(prior_dict["sqrt_corr_kernels"]*"LCO.jdat")
end

@everywhere begin
    function grab_spectrum(fname)
        return readdlm(fname)[:,2].-1
    end

    function get_star(intup)
        (fname, rv, fibernum) = intup
        rval = indInt(rv)
        tval = indTenth(rv)

        hres_spec = grab_spectrum(fname)

        # this could be more efficient, except for that we are going to 
        # batch with many fiber numbers at some point
        Ksp = if adjfibindx>300
            deserialize(prior_dict["LSF_mat_LCO"]*"$(tval)_"*lpad(adjfibindx-300,3,"0")*".jdat");
        else
            deserialize(prior_dict["LSF_mat_APO"]*"$(tval)_"*lpad(adjfibindx,3,"0")*".jdat");
        end
        nvecLSF = dropdims(sum(Ksp,dims=2),dims=2);
        lres_spec = (Ksp*hres_spec)./nvecLSF;

        return ShiftedArrays.circshift(lres_spec,rval)
    end
    
    function take_draw(ovtup; skycont_only = false, no_sky = false, dibs_on=true, correlated_noise=false, caching=true, dib_center_lambda_lst=[15273], inject_cache_dir="./inject_local_cache",cache_dir="./local_cache")
        argtup, starCont_indx, skyCont_indx, starscale, fname, rv, ew, λ0, sigma, injectindx, injectfiber = ovtup
        ival = argtup[1]
        intup = argtup[2:end]
        (release_dir,redux_ver,tele,field,plate,mjd,fiberindx) = intup

        teleind = (tele[1:6] == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiberindx
        gain = (tele[1:6] == "lco25m") ? 3.0 : 1.9

        # Get Throughput Fluxing 
        fluxingcache = cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir)
        if !isfile(fluxingcache)
            dirName = splitdir(fluxingcache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            getAndWrite_fluxing(release_dir,redux_ver,tele,field,plate,mjd,cache_dir=cache_dir)
        end
        
        skycacheSpec = cache_skynameSpec(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
        if (isfile(skycacheSpec) & caching)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skycacheSpec)
            starscalesky,framecnts,a_relFlux,b_relFlux,c_relFlux,cartVisit,ingest_bit = metaexport
        else
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
            starscalesky,framecnts,a_relFlux,b_relFlux,c_relFlux,cartVisit,ingest_bit = metaexport
            if caching
                dirName = splitdir(skycacheSpec)[1]
                if !ispath(dirName)
                    mkpath(dirName)
                end
                serialize(skycacheSpec,[fvec, fvarvec, cntvec, chipmidtimes, metaexport])
            end
        end
        outup = (release_dir,redux_ver,tele*lpad(injectindx,7,"0")*"i",field,plate,mjd,injectfiber)
        starcache = cache_starname(tele*lpad(injectindx,7,"0")*"i",field,plate,mjd,injectfiber,cache_dir=cache_dir,inject_cache_dir=inject_cache_dir)
        
        if !isfile(starcache)
            starcomp = starContSamples[:,starCont_indx]
            medstar = nanzeromedian(starcomp)
            starcomp .*= (starscale/medstar)
            starLines = get_star((fname, rv, adjfibindx))
            starcomp .*= (1 .+ starLines)
            if dibs_on
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
            end

            # add only sky cont
            if skycont_only
                starcomp .+= view(skyContSamples,:,skyCont_indx)
            end

            # seed the Poisson draw
            rng = MersenneTwister(mjd+adjfibindx) # seed on mjd and fiber index so unique
            if any(isnan.(starcomp))
                outfvec = NaN*ones(length(starcomp))
            elseif (no_sky | skycont_only)
                outfvec = pois_rand.(rng,starcomp*gain)/gain
            elseif correlated_noise
                sqrt_cov_ref = sqrtMats[tele[1:3]]
                renorm = Diagonal(sqrt_nan.(starcomp./gain))
                outfvec = fvec .+ starcomp .+ (renorm*sqrt_cov_ref)*randn(rng,length(fvec))
            else
                outfvec = fvec .+ pois_rand.(rng,starcomp*gain)/gain
            end

            outcntvec = cntvec
            outcntvec[.!chebmsk_exp] .= 0

            framecnts = maximum(outcntvec)
            simplemsk = (outcntvec.==framecnts)
            
            starscale_o = if count(simplemsk .& (.!isnan.(outfvec)))==0
                NaN
            else
                abs(nanzeromedian(outfvec[simplemsk]))
            end

            # add variance of sky observations to variance of poisson draw of the models
            outvar = if (no_sky | skycont_only)
                zeros(length(fvarvec))
            else
                fvarvec
            end
            outvar .+= starcomp./gain # add poisson noise of the model

            # just using chipmidtimes from sky exposure, won't need to use in inject WU
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

## Read in Synthetic Stellar Models
flist = sort(glob("*/*.dat",prior_dict["korg_run_path"]));
indxpassing = map(x->parse(Int,split(split(x,"_")[end],".")[1]),flist);

param_grid_file = prior_dict["korg_run_path"]*"korg_grid_params.h5"
Teffi = h5read(param_grid_file,"Teff")
loggi = h5read(param_grid_file,"logg")
m_hi = h5read(param_grid_file,"m_h")
vmici = h5read(param_grid_file,"vmic")
a_mi = h5read(param_grid_file,"a_m")
c_mi = h5read(param_grid_file,"c_m")
abundancesi = h5read(param_grid_file,"abundances")
converged_flagi = h5read(param_grid_file,"converged_flag");

Teff = Teffi[indxpassing,..]
logg = loggi[indxpassing,..]
m_h = m_hi[indxpassing,..]
vmic = vmici[indxpassing,..]
a_m = a_mi[indxpassing,..]
c_m = c_mi[indxpassing,..]
abundances = abundancesi[indxpassing,..]
converged_flag = converged_flagi[indxpassing,..]

## Read in Relevant Sky Observations (in adjfibindx)
ntuplst = deserialize(prior_dict["sky_runlist"]*lpad(adjfibindx,3,"0")*".jdat");

## Generate Injection Parameters
rng = MersenneTwister(rnd_seed)
injectfiber = fiberindx*ones(Int,nsamp) # this should not be adjusted
injectindx = 1:nsamp
sky_tup = StatsBase.sample(rng,ntuplst,nsamp)
starCont_indx = rand(rng,1:size(starContSamples,2),nsamp);
skyCont_indx = rand(rng,1:size(skyContSamples,2),nsamp); # only matters if skycont_only is true
starscale_lst = StatsBase.sample(rng,2 .*starscale_red[starscale_red.>0],nsamp)
RV = rand(rng,Uniform(RV_range_pix...),nsamp);
korgindx = rand(rng,1:length(flist),nsamp);
fname_list = flist[korgindx]
## need to add second DIB to a loop
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

itarg = Iterators.zip(sky_tup,starCont_indx,skyCont_indx,starscale_lst,fname_list,RV,eachrow(hcat(dib_ew...)),eachrow(hcat(dib_lam...)),eachrow(hcat(dib_sig...)),injectindx,injectfiber);

## Save Injection Parameters to Disk
fname = prior_dict["out_dir"]*"inject_params.h5"
dirName = splitdir(fname)[1]
if !ispath(dirName)
    mkpath(dirName)
end
f = h5open(fname,"w")
write(f,"injectfiber",injectfiber)
write(f,"injectindx",collect(injectindx))
write(f,"starscale_lst",starscale_lst)
write(f,"RV",RV)
write(f,"korgindx",korgindx)
write(f,"fname_list",fname_list)
write(f,"Teff",Teff[korgindx,..])
write(f,"logg",logg[korgindx,..])
write(f,"m_h",m_h[korgindx,..])
write(f,"vmic",vmic[korgindx,..])
write(f,"a_m",a_m[korgindx,..])
write(f,"c_m",c_m[korgindx,..])
write(f,"abundances",abundances[korgindx,..])
write(f,"converged_flag",converged_flag[korgindx,..])
for (dibind, dib_center_lambda) in enumerate(dib_center_lambda_lst)
    write(f,"dib_sig_$(dib_center_lambda)",dib_sig[dibind])
    write(f,"dib_lam_$(dib_center_lambda)",dib_lam[dibind])
    write(f,"dib_ew_$(dib_center_lambda)",dib_ew[dibind])
end
close(f)

serialize(prior_dict["out_dir"]*"inject_sky_tup.jdat",sky_tup);
serialize(prior_dict["out_dir"]*"inject_skyCont_indx.jdat",skyCont_indx);

println("All Star Models Used Converged: ", sum(converged_flag[korgindx,..]) == length(converged_flag[korgindx,..]))

## Create an injection test series
@everywhere take_draw_partial(ovtup) = take_draw(ovtup,skycont_only=skycont_only,no_sky=no_sky,dibs_on=dibs_on,correlated_noise=correlated_noise,dib_center_lambda_lst=dib_center_lambda_lst,inject_cache_dir=prior_dict["inject_cache_dir"],cache_dir=prior_dict["local_cache"])
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
