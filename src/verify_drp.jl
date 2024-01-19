##### Verify APOGEE DRP Products from Release are Production Ready

import Pkg
using InteractiveUtils; versioninfo()
Pkg.activate("../"); Pkg.instantiate(); Pkg.precompile()

using FITSIO, StatsBase, ProgressMeter, Distributed, Serialization, Glob, DelimitedFiles

src_dir = abspath("./")
include(src_dir*"/fileNameHandling.jl")

outdir = "../../../2024_01_09/outlists/"
if !ispath(outdir)
    mkpath(outdir)
end

release_dir = "sdsswork/mwm"
redux_ver = "1.2"

release_dir_n = replace(release_dir,"/"=>"_")
redux_ver_n = replace(redux_ver,"."=>"_")

check_ap1d = false
check_apCframes = false
check_exp = false
check_flux = true

### Ingest allVisit File
tele_try_list = ["apo25m","lco25m"]
tele_list = []
for tele in tele_try_list
    allVisitpath = getUtahBase(release_dir, redux_ver)*"summary/allVisit-$(redux_ver)-$(tele).fits"
    # println(allVisitpath)
    if isfile(allVisitpath)
        push!(tele_list,tele)
    end
end

if length(tele_list) == 0
    println("No Summary Files for Either Instrument?")
end

TELESCOPE_LST = []
FIELD_LST = []
PLATE_LST = []
MJD_LST = []
FIBERID_LST = []
for tele in tele_list
    f = FITS(getUtahBase(release_dir, redux_ver)*"summary/allVisit-$(redux_ver)-$(tele).fits")
    push!(TELESCOPE_LST,read(f[2],"TELESCOPE"))
    push!(FIELD_LST,read(f[2],"FIELD"))
    push!(PLATE_LST,read(f[2],"PLATE"))
    push!(MJD_LST,read(f[2],"MJD"))
    push!(FIBERID_LST,read(f[2],"FIBERID"))
    close(f)
end

TELESCOPE = vcat(TELESCOPE_LST...)
FIELD = vcat(FIELD_LST...)
PLATE = vcat(PLATE_LST...)
MJD = vcat(MJD_LST...)
FIBERID = vcat(FIBERID_LST...)

### Boot Up Workers
addprocs(32)
@everywhere begin
    import Pkg
    Pkg.activate("../")

    using FITSIO, StatsBase, ProgressMeter, Distributed, Serialization, Glob, DelimitedFiles

    src_dir = abspath("./")
    include(src_dir*"/fileNameHandling.jl")
end

println("##########################################################################################")
println("############## Checking APOGEE DRP $(release_dir), version $(redux_ver) products ####################")
println("##########################################################################################")

### Get Stars to Run (write apMADGICS.jl run files)
### Check Last Entry of ap1D Files

@everywhere begin
    function ap1D_check(intup)
        runnum, release_dir,redux_ver,tele,field,plate,mjd,fiberindx=intup

        plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
        frame_lst = getFramesFromPlate(plateFile)

        outlst = []
        allgood = true
        for imid in frame_lst
            for chip in ["a","b","c"]
                fname = build_framepath(release_dir,redux_ver,tele,mjd,imid,chip)
                try 
                    f = FITS(fname)
                    temp = read(f[5],:,fiberindx)
                    close(f)
                    allgood &= true
                catch
                    allgood &= false
                end
            end
        end
        push!(outlst,allgood)
        return outlst
    end

    function get_bad_ap1D(intup)
        runnum, release_dir,redux_ver,tele,field,plate,mjd,fiberindx=intup

        plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
        frame_lst = getFramesFromPlate(plateFile)

        outlst = []
        allgood = true
        for imid in frame_lst
            for chip in ["a","b","c"]
                fname = build_framepath(release_dir,redux_ver,tele,mjd,imid,chip)
                try 
                    f = FITS(fname)
                    temp = read(f[5],:,fiberindx)
                    close(f)
                    allgood &= true
                catch
                    allgood &= false
                    push!(outlst,fname)
                end
            end
        end
        return outlst
    end
end

if check_ap1d
    println("##### Checking the ap1D Files #####")

    starlst = []
    badlst = []
    for telematch in tele_list
        @showprogress for fiber in 1:300
            teleind = (telematch == "lco25m") ? 2 : 1
            adjfibindx = (teleind-1)*300 + fiber
            mskfib = (FIBERID.==(301-fiber)) .& (TELESCOPE .== telematch);
            nspectra = count(mskfib)
            # println((telematch, fiber, nspectra))
            if nspectra > 1
                subiter = Iterators.zip(1:nspectra,repeat([release_dir],nspectra),repeat([redux_ver],nspectra),TELESCOPE[mskfib],FIELD[mskfib],string.(parse.(Int,PLATE[mskfib])),MJD[mskfib],fiber*ones(Int,nspectra))
                tout = pmap(ap1D_check,subiter)
                outcheck = convert(Vector{Bool},vcat(tout...))
                serialize(outdir*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat",collect(subiter)[outcheck])
                if count(.!outcheck)>0
                    serialize(outdir*"$(release_dir_n)_$(redux_ver_n)_star_badap1D_lst_"*lpad(adjfibindx,3,"0")*".jdat",collect(subiter)[.!outcheck])
                end
                push!(starlst,count(outcheck))
            end
            flush(stdout)
        end
    end

    badlst = []
    for telematch in tele_list
        @showprogress for fiber in 1:300
            teleind = (telematch == "lco25m") ? 2 : 1
            adjfibindx = (teleind-1)*300 + fiber
            fname = outdir*"star_bad_lst_"*lpad(adjfibindx,3,"0")*".jdat"
            if isfile(fname)
                subiter = deserialize(fname)
                badframes = get_bad_ap1D.(subiter)
                push!(badlst,badframes)
            end
        end
    end

    badfilevec = unique(vcat(map(x->vcat(x...),badlst)...))

    if length(badfilevec)>0
        println("Total bad ap1D files: ",length(badfilevec))
        writedlm(outdir*"$(release_dir_n)_$(redux_ver_n)_bad_ap1D.txt",badfilevec,',')
    end
end

### Check apCframe Files (change from exist to last entry)

if check_apCframes
    println("##### Checking the apCframe Files #####")

    apCframes = sort(glob("*/*/*/*/apCframe-*-*.fits",getUtahBase(release_dir, redux_ver)*"visit"))

    @everywhere begin
        function check_file_size(fname)
            return stat(fname).size
        end

        function apCframe_check(fname)
            allgood = true
            try 
                f = FITS(fname)
                close(f)
                allgood &= true
            catch
                allgood &= false
            end
            return allgood
        end
    end

    fcheck = check_file_size.(apCframes)
    println("Total apCframe files that are Zero bytes: ",count(fcheck.==0))

    # fig = plt.figure(figsize=(6,6),dpi=150)
    # ax = fig.add_subplot(1,1,1)
    # ax.hist(fcheck,bins=100);
    # ax.set_yscale("log")
    # fig.savefig(outdir*"$(release_dir)_$(redux_ver)_apCframe_size.png", bbox_inches="tight", pad_inches=0.1);
    # plt.close()

    scheck = @showprogress pmap(apCframe_check,apCframes);
    ocheck = vcat(scheck...)

    if count(.!ocheck)>0
        println("Total bad apCframe files: ",count(.!ocheck))
        writedlm(outdir*"$(release_dir_n)_$(redux_ver_n)_bad_apCframe.txt",fnames[.!ocheck],',')
    end
end

@everywhere begin
    function exp_check(intup)
        runnum,release_dir,redux_ver,tele,field,plate,mjd,fiberindx=intup

        plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
        frame_lst = getFramesFromPlate(plateFile)
        firstVisitExp =  minimum(parse.(Int,frame_lst))
        expFile = build_expPath(release_dir,redux_ver,tele,mjd)
        if !isfile(expFile)
            return 3
        end
        try
            f = FITS(expFile)
            MED = read(f[2],"MED")
            close(f)
            return 0
        catch
            return 1
        end
    end
end

run_lsts = []
run_lens = []
for telematch in tele_list
    for fiber in 1:300
        teleind = (telematch == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiber
        run_per_fiber = deserialize(outdir*"star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat")
        push!(run_lsts,run_per_fiber)
        push!(run_lens,length(run_per_fiber))
    end
end
run_lst = vcat(run_lsts...)

if check_exp
    println("##### Checking the MJDexp Files #####")
    pout = @showprogress pmap(exp_check,run_lst)

    if count((pout.!=0))>0
        mjd_bad = unique(map(x->x[7],run_lst[pout.!=0]))
        println("Missing exp files: ",length(mjd_bad))
        println("Uknown errors for exp files: ",count(pout.!=0)-count((pout.==3)))
        writedlm(outdir*"$(release_dir_n)_$(redux_ver_n)_bad_exp.txt",mjd_bad,',')
    end
end

@everywhere begin
    function test_apFlux(intup;ret_file=false)
        runnum,release_dir,redux_ver,tele,field,plate,mjd,fiberindx=intup
        
        plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
        frame_lst = getFramesFromPlate(plateFile)
        firstVisitExp =  minimum(parse.(Int,frame_lst))
        expFile = build_expPath(release_dir,redux_ver,tele,mjd)
        
        if !isfile(expFile)
            return 3
        end
        f = FITS(expFile)
        EXPNUM = read(f[2],"NUM")
        IMAGETYP = read(f[2],"IMAGETYP")
        CARTID = read(f[2],"CARTID")
        close(f)
        
        firstVisitInd = findfirst(EXPNUM.==firstVisitExp)
        cartVisit = CARTID[firstVisitInd]

        expIndex_before = findlast((IMAGETYP.=="DomeFlat") .& (EXPNUM .< firstVisitExp))
        expIndex_after = findfirst((IMAGETYP.=="DomeFlat") .& (EXPNUM .> firstVisitExp))
        
        # This is great logic, if the DRP did not only do the first two in the FPI era
        # I am leaving it and overriding it because I want to be able to revert to this
        expIndex = if isnothing(expIndex_before) & isnothing(expIndex_after)
            -1
        elseif !isnothing(expIndex_before) && (CARTID[expIndex_before] == cartVisit)
            expIndex_before
        elseif !isnothing(expIndex_after) && (CARTID[expIndex_after] == cartVisit)
            expIndex_after
        else
            -2
        end
        
        for off in vcat(0,-1:-1:-3,1:3)
            if (1 <= expIndex+off <= length(CARTID)) && isfile(build_apFluxPath(release_dir,redux_ver,tele,mjd,"a",EXPNUM[expIndex+off])) && (CARTID[expIndex+off] == cartVisit)
                expIndex+=off
                break
            end
        end
        
        # This is to handle bad DRP behavior upstream
        if all(CARTID.==0)
            allDomes = findall(IMAGETYP.=="DomeFlat")
            testinds = allDomes[(allDomes .<= expIndex)]
            for testind in testinds[end:-1:1]
                fluxFile = build_apFluxPath(release_dir,redux_ver,tele,mjd,"a",EXPNUM[testind])
                if isfile(fluxFile)
                    expIndex = testind
                    break
                else
                    expIndex = -3
                end
            end
        end
            
        if expIndex < 0
            return expIndex
        end
        expectDome = EXPNUM[expIndex]

        allgood = 0
        for chip in ["a","b","c"]
            fluxFile = build_apFluxPath(release_dir,redux_ver,tele,mjd,chip,expectDome)
            if !isfile(fluxFile)
                if ret_file
                    return fluxFile
                else
                    return 4
                end
            end
            try
                f = FITS(fluxFile)
                thrpt = read(f[3])
                lastExt = read(f[5])
                close(f)
                allgood = 0
            catch
                allgood = 1
            end
        end
        return allgood
    end

    g(x) = test_apFlux(x,ret_file=true)

    function divide_vector(input_vector, lengths)
        if sum(lengths) != length(input_vector)
            throw(ArgumentError("Sum of lengths must be equal to the length of the input vector"))
        end
        divided_vectors = []
        start_indices = [0; cumsum(lengths)]
        for i in 1:length(lengths)
            push!(divided_vectors,input_vector[start_indices[i]+1:start_indices[i+1]])
        end
        return divided_vectors
    end
end

if check_flux
    println("##### Checking the apFlux Files #####")
    pout = @showprogress pmap(test_apFlux,run_lst)
    pout1 = @showprogress pmap(g,run_lst);

    if count((pout.!=0))>0
        println("Visits lost to missing exp files: $(count(pout.==3)) , $(100*count(pout.==3)/length(pout))%")
        println("Visits lost because no flat taken at all: $(count(pout.==-1)) , $(100*count(pout.==-1)/length(pout))%")
        println("Visits lost because no flat taken with CARTID: $(count(pout.==-2)) , $(100*count(pout.==-2)/length(pout))%")
        println("Visits lost to missing flux files: $(count((pout.==4).|(pout.==-3))) , $(100*count((pout.==4).|(pout.==-3))/length(pout))%")
        println("Unknown errors for flux files: ",count(pout.!=0)-count((pout.==3) .| (pout.==4) .| (pout.==-1) .| (pout.==-2) .| (pout.==-3)))

        writedlm(outdir*"$(release_dir_n)_$(redux_ver_n)_missing_apFlux.txt",pout1[(pout.==4).|(pout.==-3)],',')
        
        missFluxing = divide_vector(pout.!=0, run_lens)
        global ind_loc = 1
        for telematch in tele_list
            for fiber in 1:300
                teleind = (telematch == "lco25m") ? 2 : 1
                adjfibindx = (teleind-1)*300 + fiber
                serialize(outdir*"$(release_dir_n)_$(redux_ver_n)_star_msk_fluxing_lst_"*lpad(adjfibindx,3,"0")*".jdat",missFluxing[ind_loc])
                global ind_loc+=1
            end
        end

        mjd_bad = unique(map(x->x[7],run_lst[(pout.==4).|(pout.==-3)]))
        println("MJDs with missing apFlux Files: ",length(mjd_bad))
        writedlm(outdir*"$(release_dir_n)_$(redux_ver_n)_bad_mjd_missing_apFlux.txt",mjd_bad,',')
    end
end

# apply  any other masking to the star list and generate final star lists
for telematch in tele_list
    for fiber in 1:300
        teleind = (telematch == "lco25m") ? 2 : 1
        adjfibindx = (teleind-1)*300 + fiber
        star_input = deserialize(outdir*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat")
        msk = deserialize(outdir*"$(release_dir_n)_$(redux_ver_n)_star_msk_fluxing_lst_"*lpad(adjfibindx,3,"0")*".jdat")
        serialize(outdir*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_msked"*lpad(adjfibindx,3,"0")*".jdat",star_input[.!msk])
    end
end

rmprocs(workers())

