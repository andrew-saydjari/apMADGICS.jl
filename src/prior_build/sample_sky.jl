## This is script grabs a bunch of sky fiber spectra and decomposes them into continuum and line components, to serve as samples for building the sky prior
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
    using SortFilters, BasisFunctions, Random
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    # runlist_range = 295 #1:600 #295, 245, 335, 101

    # Prior Dictionary
    prior_dict = Dict{String,String}()

    # Input List (not really a prior, but an input file we search for stars conditioned on)
    prior_dict["sky_runlist"] = prior_dir*"2024_02_21/outlists/sky/dr17_dr17_sky_input_lst_plate_msked_"

    # Data for Detector Cals (not really a prior, but an input the results depend on in detail)
    prior_dict["chip_fluxdep_err_correction"] = src_dir*"data/chip_fluxdep_err_correction.jdat"
    prior_dict["medframes_APO"] = src_dir*"data/medframes_apo.jdat" # last made 2023_04_01 by AKS
    prior_dict["medframes_LCO"] = src_dir*"data/medframes_lco.jdat" # last made 2023_04_07 by AKS
end

@everywhere begin
    # I should revisit the error bars in the context of chi2 versus frame number trends
    global err_correct_Dict = deserialize(prior_dict["chip_fluxdep_err_correction"])

    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg)
    
    c = 299792.458; # in km/s
    delLog = 6e-6; 
    # pixscale = (10^(delLog)-1)*c; # only a linear approx, use more accurate formula

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
    function ingest_skyVisit_stack(argtup;cache_dir="./local_cache_sky")
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

        skyLineCache = cache_skynameSpec(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
        if !isfile(skyLineCache)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
            dirName = splitdir(skyLineCache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(skyLineCache,[fvec, fvarvec, cntvec, chipmidtimes, metaexport])
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
end
@everywhere begin
    function sky_smooth_wrapper(argtup,chebmsk_exp,Vpoly_scaled;telluric_div=false,cache_dir="./local_cache_sky")
        try
            release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

            skyLineCache_tellDiv = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvecN = deserialize(skyLineCache_tellDiv)
            simplemsk = (cntvec.==maximum(cntvec)) .& chebmsk_exp;
            if telluric_div
                fnew = sky_smooth_fit(fvec,fvarvec,simplemsk,Vpoly_scaled)
            else
                fnew = sky_smooth_fit(fvec,fvarvec,simplemsk,Vpoly_scaled).*telvecN
            end
            return fnew
        catch
            release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

            skyLineCache_tellDiv = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
            println("Failed to load ",skyLineCache_tellDiv); flush(stdout)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvecN = deserialize(skyLineCache_tellDiv)
            simplemsk = (cntvec.==maximum(cntvec)) .& chebmsk_exp;
            if telluric_div
                fnew = sky_smooth_fit(fvec,fvarvec,simplemsk,Vpoly_scaled)
            else
                fnew = sky_smooth_fit(fvec,fvarvec,simplemsk,Vpoly_scaled).*telvecN
            end
            return fnew
        end
    end

    function sky_line_wrapper(argtup,skycont;telluric_div=false,cache_dir="./local_cache_sky")
        ival = argtup[1]
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=telluric_div,cache_dir=cache_dir)
        fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skyLineCache) #ignore telvecN since not using even in tell_div case
        simplemsk = (cntvec.==maximum(cntvec));
        
        fnew = fvec.-skycont[:,ival]
        fnew[.!simplemsk].=0
        return fnew
    end

    function sky_msk_wrapper(argtup;cache_dir="./local_cache_sky")
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache = cache_skynameSpec(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
        fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skyLineCache)
        simplemsk = (cntvec.==maximum(cntvec));
        return simplemsk, fvarvec
    end

    function sky_smooth_fit(outvec,outvar,simplemsk,Vpoly_scaled)
        ## Select data for use (might want to handle mean more generally)
        wave_obs = wavetarg[simplemsk]
        Xd_obs0 = outvec[simplemsk];

        tmsk = ret_qlines(Xd_obs0,wave_obs)
        Xd_obs = outvec[simplemsk][tmsk];

        ## Set up residuals prior
        A = Diagonal(outvar[simplemsk][tmsk]);
        Ainv = Diagonal(1 ./outvar[simplemsk][tmsk]);

        ## Set up priors
        V_smooth_c = Vpoly_scaled
        V_smooth_r = V_smooth_c[simplemsk,:][tmsk,:]

        # Compute sky line/continuum separation
        Vcomb = V_smooth_r
        Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
        x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, (V_smooth_r, ), (V_smooth_c, ))

        return x_comp_lst[1]
    end
end

@everywhere begin
    # contScale is a tuning parameter we might want to investigate
    function get_sky_samples(adjfibindx;contscale=5e2,loc_parallel=false,seed=2023)
        
        ntuplst = deserialize(prior_dict["sky_runlist"]*lpad(adjfibindx,3,"0")*".jdat")

        medframes = if adjfibindx>300
            deserialize(prior_dict["medframes_LCO"]);
        else
            deserialize(prior_dict["medframes_APO"]);
        end

        Vpoly, msknall = generate_poly_prior(adjfibindx,medframes)
        Vpoly_scaled = contscale*Vpoly
        chebmsk = .!msknall
        chebmsk_exp = expand_msk(chebmsk;rad=12);

        ### Stack and Strip Sky From Visits
        pout = if loc_parallel
            @showprogress pmap(ingest_skyVisit_stack,ntuplst);
        else
            map(ingest_skyVisit_stack,ntuplst);
        end

        ### Save samples with usual telluric contributions included
        savename = "sky_prior_disk/skycont_"*lpad(adjfibindx,3,"0")*".jdat"
        sky_smooth_wrapper_bound(argtup) = sky_smooth_wrapper(argtup,chebmsk_exp,Vpoly_scaled)
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_smooth_wrapper_bound,ntuplst);
            else
                map(sky_smooth_wrapper_bound,ntuplst);
            end
            skycont = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skycont[:,i].=pout[i]
            end
            
            dirName = splitdir(savename)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(savename,skycont)
        else
            global skycont = deserialize(savename)
        end

        savename = "sky_prior_disk/skyline_"*lpad(adjfibindx,3,"0")*".jdat"
        sky_line_wrapper_bound(argtup) = sky_line_wrapper(argtup,skycont)
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_line_wrapper_bound,ntuplst);
            else
                map(sky_line_wrapper_bound,ntuplst);
            end
            skyline = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skyline[:,i].=pout[i]
            end
            serialize(savename,skyline)
        end

        savename = "sky_prior_disk/skymsk_"*lpad(adjfibindx,3,"0")*".jdat"
        savename1 = "sky_prior_disk/skyvar_"*lpad(adjfibindx,3,"0")*".jdat"
        if !(isfile(savename) & isfile(savename1))
            pout = if loc_parallel 
                @showprogress pmap(sky_msk_wrapper,ntuplst);
            else
                map(sky_msk_wrapper,ntuplst);
            end
            skymsk = zeros(8700,size(pout,1));
            skyvar = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skymsk[:,i].=pout[i][1]
                skyvar[:,i].=pout[i][2]
            end
            serialize(savename,skymsk)
            serialize(savename1,skyvar)
        end

        ### Save samples of tell-free sky decomposition for building Tfun/starCont prior
        sky_smooth_wrapper_tell(argtup) = sky_smooth_wrapper(argtup,chebmsk_exp,Vpoly_scaled,telluric_div=true)
        sky_line_wrapper_tell(argtup) = sky_line_wrapper(argtup,skycont,telluric_div=true)

        savename = "sky_prior_disk/skycont_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_smooth_wrapper_tell,ntuplst);
            else
                map(sky_smooth_wrapper_tell,ntuplst);
            end
            skycont = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skycont[:,i].=pout[i]
            end
            
            dirName = splitdir(savename)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(savename,skycont)
        end

        savename = "sky_prior_disk/skyline_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_line_wrapper_tell,ntuplst);
            else
                map(sky_line_wrapper_tell,ntuplst);
            end
            skyline = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skyline[:,i].=pout[i]
            end
            serialize(savename,skyline)
        end

        # this is identical... just saving under a new name, but it is cheap
        savename = "sky_prior_disk/skymsk_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        savename1 = "sky_prior_disk/skyvar_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        if !(isfile(savename) & isfile(savename1))
            pout = if loc_parallel 
                @showprogress pmap(sky_msk_wrapper,ntuplst);
            else
                map(sky_msk_wrapper,ntuplst);
            end
            skymsk = zeros(8700,size(pout,1));
            skyvar = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skymsk[:,i].=pout[i][1]
                skyvar[:,i].=pout[i][2]
            end
            serialize(savename,skymsk)
            serialize(savename1,skyvar)
        end
        serialize("sky_prior_disk/chebmsk_exp_"*lpad(adjfibindx,3,"0")*".jdat",chebmsk_exp)
    end
end

# get_sky_samples(runlist_range,loc_parallel=true) #295, 335, 450, 101
@showprogress pmap(get_sky_samples,1:600) # 13ish hours on 4 np nodes (this included the svd... which I have moved out now)
