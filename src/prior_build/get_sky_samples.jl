## This is script grabs a bunch of sky fiber spectra and decomposes them into continuum and line components, to serve as samples for building the sky prior
# Author - Andrew Saydjari, CfA
import Pkg; using Dates; t0 = now(); t_then = t0;
using InteractiveUtils; versioninfo()
Pkg.activate("./"); Pkg.instantiate(); Pkg.precompile()
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Package activation took $dt"); t_then = t_now; flush(stdout)
using Distributed, SlurmClusterManager, Suppressor, DataFrames
addprocs(SlurmManager())
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker allocation took $dt"); t_then = t_now; flush(stdout)
activateout = @capture_out begin
    @everywhere begin
        import Pkg
        Pkg.activate("./")
    end
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker activation took $dt"); t_then = t_now; flush(stdout)
@everywhere begin
    using FITSIO, Serialization, HDF5, LowRankOps, EllipsisNotation, ShiftedArrays
    using Interpolations, SparseArrays, ParallelDataTransfer, AstroTime, Suppressor
    using ThreadPinning

    prior_dir = "/uufs/chpc.utah.edu/common/home/u6039752/scratch1/working/"
    src_dir = "../../apMADGICS.jl/"
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
    
    using StatsBase, LinearAlgebra, ProgressMeter
    using BLISBLAS
    BLAS.set_num_threads(1)
end
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("Worker loading took $dt"); t_then = t_now; flush(stdout)

# Task-Affinity CPU Locking in multinode SlurmContext
slurm_cpu_lock()
println(BLAS.get_config()); flush(stdout)
t_now = now(); dt = Dates.canonicalize(Dates.CompoundPeriod(t_now-t_then)); println("CPU locking took $dt"); t_then = t_now; flush(stdout)

using LibGit2; git_branch, git_commit = initalize_git(src_dir); @passobj 1 workers() git_branch; @passobj 1 workers() git_commit

@everywhere begin
    wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125) #first argument is start, revert fix to enable 1.6 compat
    minw, maxw = extrema(wavetarg)
    
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
    function ingest_skyVisit_stack(argtup;cache_dir="./local_cache_sky")
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache = cache_skynameSpec(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
        if !isfile(skyLineCache)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx)
            dirName = splitdir(skyLineCache)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(skyLineCache,[fvec, fvarvec, cntvec, chipmidtimes, metaexport])
        end

        skyLineCache_tellDiv = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
        if !isfile(skyLineCache_tellDiv)
            fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvec = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx,telluric_div=true)
            dirName = splitdir(skyLineCache_tellDiv)[1]
            if !ispath(dirName)
                mkpath(dirName)
            end
            serialize(skyLineCache_tellDiv,[fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvec])
        end
    end
end
@everywhere begin
    function sky_smooth_wrapper(argtup;telluric_div=false,cache_dir="./local_cache_sky")
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache_tellDiv = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=true,cache_dir=cache_dir)
        fvec, fvarvec, cntvec, chipmidtimes, metaexport, telvecN = deserialize(skyLineCache_tellDiv)
        simplemsk = (cntvec.==maximum(cntvec)) .& chebmsk_exp;
        if telluric_div
            fnew = sky_smooth_fit(fvec,fvarvec,simplemsk)
        else
            fnew = sky_smooth_fit(fvec,fvarvec,simplemsk).*telvecN
        end
        return fnew
    end

    function sky_line_wrapper(argtup;telluric_div=false,cache_dir="./local_cache_sky")
        ival = argtup[1]
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache = cache_skynameSpec(tele,field,plate,mjd,fiberindx,telluric_div=telluric_div,cache_dir=cache_dir)
        fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skyLineCache) #ignore telvecN since not using even in tell_div case
        simplemsk = (cntvec.==maximum(cntvec));
        
        fnew = fvec.-skycont[:,ival] #### UHHHHH IS THIS A GROSS GLOBAL???? YES :( FIXME
        fnew[.!simplemsk].=0
        return fnew
    end

    function sky_msk_wrapper(argtup;cache_dir="./local_cache_sky")
        release_dir, redux_ver, tele, field, plate, mjd, fiberindx = argtup[2:end]

        skyLineCache = cache_skynameSpec(tele,field,plate,mjd,fiberindx,cache_dir=cache_dir)
        fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skyLineCache)
        simplemsk = (cntvec.==maximum(cntvec));
        return simplemsk
    end

    function sky_smooth_fit(outvec,outvar,simplemsk)
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
    function get_sky_samples(adjfibindx,contscale=5e2,loc_parallel=false)
        
        ntuplst = deserialize("../2023_05_21/sky_input_lists/sky_input_lst_"*lpad(adjfibindx,3,"0")*".jdat")

        if adjfibindx>300
            global medframes = deserialize(prior_dir*"2023_04_07/medframes_lco.jdat");
        else
            global medframes = deserialize(prior_dir*"2023_04_01/medframes.jdat");
        end

        global Vpoly, msknall = generate_poly_prior(adjfibindx)
        global Vpoly_scaled = contscale*Vpoly
        global chebmsk = .!msknall
        global chebmsk_exp = expand_msk(chebmsk;rad=12);

        ### Stack and Strip Sky From Visits
        pout = if loc_parallel
            @showprogress pmap(ingest_skyVisit_stack,ntuplst);
        else
            map(ingest_skyVisit_stack,ntuplst);
        end

        ### Save samples with usual telluric contributions included
        savename = "sky_prior_disk/skycont_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_smooth_wrapper,ntuplst);
            else
                map(sky_smooth_wrapper,ntuplst);
            end
            global skycont = zeros(8700,size(pout,1));
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
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_line_wrapper,ntuplst);
            else
                map(sky_line_wrapper,ntuplst);
            end
            global skyline = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skyline[:,i].=pout[i]
            end
            serialize(savename,skyline)
        else
            global skyline = deserialize(savename)
        end

        savename = "sky_prior_disk/skymsk_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel 
                @showprogress pmap(sky_msk_wrapper,ntuplst);
            else
                map(sky_msk_wrapper,ntuplst);
            end
            global skymsk = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skymsk[:,i].=pout[i]
            end
            serialize(savename,skymsk)
        else
            global skymsk = deserialize(savename)
        end

        ### Save samples of tell-free sky decomposition for building Tfun/starCont prior
        sky_smooth_wrapper_tell(argtup) = sky_smooth_wrapper(argtup,telluric_div=true)
        sky_line_wrapper_tell(argtup) = sky_line_wrapper(argtup,telluric_div=true)

        savename = "sky_prior_disk/skycont_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_smooth_wrapper_tell,ntuplst);
            else
                map(sky_smooth_wrapper_tell,ntuplst);
            end
            global skycont = zeros(8700,size(pout,1));
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

        savename = "sky_prior_disk/skyline_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel
                @showprogress pmap(sky_line_wrapper_tell,ntuplst);
            else
                map(sky_line_wrapper_tell,ntuplst);
            end
            global skyline = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skyline[:,i].=pout[i]
            end
            serialize(savename,skyline)
        else
            global skyline = deserialize(savename)
        end

        # this is identical... just saving under a new name, but it is cheap
        savename = "sky_prior_disk/skymsk_tellDiv_"*lpad(adjfibindx,3,"0")*".jdat"
        if !isfile(savename)
            pout = if loc_parallel 
                @showprogress pmap(sky_msk_wrapper,ntuplst);
            else
                map(sky_msk_wrapper,ntuplst);
            end
            global skymsk = zeros(8700,size(pout,1));
            for i=1:size(pout,1)
                skymsk[:,i].=pout[i]
            end
            serialize(savename,skymsk)
        else
            global skymsk = deserialize(savename)
        end
    end
end

get_sky_samples(295,loc_parallel=true)
# @showprogress pmap(get_sky_samples,1:600) # 13ish hours on 4 np nodes (this included the svd... which I have moved out now)
# can come back to the question of doing the svd and prior building here while the sky samples are in memory

