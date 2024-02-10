##### Verify APOGEE DRP Products from Release are Production Ready

println("##########################################################################################")
println("############################## Starting verify_drp.jl ####################################")
println("##########################################################################################")

import Pkg
using InteractiveUtils; versioninfo()
Pkg.activate("../"); Pkg.instantiate(); Pkg.precompile()

using FITSIO, StatsBase, ProgressMeter, Distributed, Serialization, Glob, DelimitedFiles, Suppressor

src_dir = abspath("./")
include(src_dir*"/fileNameHandling.jl")
include(src_dir*"/utils.jl")

outdir = "../../outlists/"
if !ispath(outdir)
    mkpath(outdir)
end
outdirsubs = [outdir*"summary/",outdir*"star/",outdir*"sky/",outdir*"tell/"]
for outsub in outdirsubs
    if !ispath(outsub)
        mkpath(outsub)
    end
end

# release_dir = "ipl-3"
# redux_ver = "1.2"

# release_dir = "dr17"
# redux_ver = "dr17"

release_dir = ARGS[1]
redux_ver = ARGS[2]

release_dir_n = replace(replace(release_dir,"/"=>"_"),"-"=>"_")
redux_ver_n = replace(redux_ver,"."=>"p")

check_ap1d = true
check_apCframes = true
check_exp = true
check_flux = true
check_plates = true
write_sky = true
write_star_plate = true
write_tell = true
check_wavecal = true
tele_try_list =  ["apo25m","lco25m"]


### Ingest allVisit File
function summary_file_by_dr(release_dir,redux_ver,dr_number)
    if dr_number==17
        return replace(getUtahBase(release_dir,redux_ver),"redux"=>"aspcap")*"synspec/allVisit-$(redux_ver)-synspec.fits"
    elseif dr_number==16
        return replace(getUtahBase(release_dir,redux_ver),"redux"=>"aspcap")*"l33/allVisit-$(redux_ver)-l33.fits"
    elseif dr_number==15
        return getUtahBase(release_dir,redux_ver)*"allVisit-l31c.2.fits"
    else
        error("No support for $(dr_number) yet, due to inhomogenous summary files")
    end
end

function ingest_allVisit_file(release_dir,redux_ver;tele_try_list=["apo25m","lco25m"])
    dr_number = if occursin("dr", release_dir)
        parse(Int, match(r"dr(\d+)", release_dir).captures[1])
    else
        -1
    end
    
    if (10 <= dr_number <= 17)
        # there is only one for the early DRs
        sum_file = summary_file_by_dr(release_dir,redux_ver,dr_number)
        
        f = FITS(sum_file)
        TELESCOPE = read(f[2],"TELESCOPE")
        FIELD = read(f[2],"FIELD")
        PLATE = read(f[2],"PLATE")
        MJD = read(f[2],"MJD")
        FIBERID = read(f[2],"FIBERID")
        close(f)
        
        msk_tele = zeros(Bool,length(TELESCOPE))
        for tele_try in tele_try_list
            msk_tele .|= (TELESCOPE .== tele_try)
        end
        
        tele_list = unique(TELESCOPE[msk_tele])

        println("Found $(length(TELESCOPE[msk_tele])) visits to run")
        
        return tele_list, TELESCOPE[msk_tele], FIELD[msk_tele], PLATE[msk_tele], MJD[msk_tele], FIBERID[msk_tele]
        ## TBD
    else
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

        println("Found $(length(TELESCOPE)) visits to run")
        
        return tele_list, TELESCOPE, FIELD, PLATE, MJD, FIBERID
    end
end

### Boot Up Workers
addprocs(32)
activateout = @capture_out begin 
    @everywhere begin
        import Pkg
        Pkg.activate("../"); Pkg.instantiate(); #comment out instantiate next time
    end
end

@everywhere begin
    using FITSIO, StatsBase, ProgressMeter, Distributed, Serialization, Glob, DelimitedFiles, ParallelDataTransfer, DataFrames

    src_dir = abspath("./")
    include(src_dir*"/fileNameHandling.jl")
    include(src_dir*"/utils.jl")
end

println("##########################################################################################")
println("############## Checking APOGEE DRP $(release_dir), version $(redux_ver) products ####################")
println("##########################################################################################")
println("##### Considering files from $(tele_try_list) #####")

tele_list, TELESCOPE, FIELD, PLATE, MJD, FIBERID = ingest_allVisit_file(release_dir,redux_ver,tele_try_list=tele_try_list)

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
                serialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat",collect(subiter)[outcheck])
                if count(.!outcheck)>0
                    serialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_badap1D_lst_"*lpad(adjfibindx,3,"0")*".jdat",collect(subiter)[.!outcheck])
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
            fname = outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_badap1D_lst_"*lpad(adjfibindx,3,"0")*".jdat"
            if isfile(fname)
                subiter = deserialize(fname)
                badframes = get_bad_ap1D.(subiter)
                push!(badlst,badframes)
            end
        end
    end

    badfilevec = unique(vcat(map(x->vcat(x...),badlst)...))
    
    println("Total bad ap1D files: ",length(badfilevec))
    if length(badfilevec)>0
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_ap1D.txt",badfilevec,',')
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
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_apCframe.txt",fnames[.!ocheck],',')
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

if check_exp
    run_lsts = []
    run_lens = []
    for telematch in tele_list
        for fiber in 1:300
            teleind = (telematch == "lco25m") ? 2 : 1
            adjfibindx = (teleind-1)*300 + fiber
            run_per_fiber = deserialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat")
            push!(run_lsts,run_per_fiber)
            push!(run_lens,length(run_per_fiber))
        end
    end
    run_lst = vcat(run_lsts...)

    println("##### Checking the MJDexp Files #####")
    pout = @showprogress pmap(exp_check,run_lst)

    if count((pout.!=0))>0
        mjd_bad = unique(map(x->x[7],run_lst[pout.!=0]))
        println("Missing exp files: ",length(mjd_bad))
        println("Uknown errors for exp files: ",count(pout.!=0)-count((pout.==3)))
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_exp.txt",mjd_bad,',')
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

        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_missing_apFlux.txt",pout1[(pout.==4).|(pout.==-3)],',')
        
        missFluxing = divide_vector(pout.!=0, run_lens)
        global ind_loc = 1
        for telematch in tele_list
            for fiber in 1:300
                teleind = (telematch == "lco25m") ? 2 : 1
                adjfibindx = (teleind-1)*300 + fiber
                serialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_msk_fluxing_lst_"*lpad(adjfibindx,3,"0")*".jdat",missFluxing[ind_loc])
                global ind_loc+=1
            end
        end

        mjd_bad = unique(map(x->x[7],run_lst[(pout.==4).|(pout.==-3)]))
        println("MJDs with missing apFlux Files: ",length(mjd_bad))
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_mjd_missing_apFlux.txt",mjd_bad,',')
    end
end

if check_ap1d | check_apCframes | check_exp | check_flux
    # apply any other masking to the star list and generate final star lists
    for telematch in tele_list
        for fiber in 1:300
            teleind = (telematch == "lco25m") ? 2 : 1
            adjfibindx = (teleind-1)*300 + fiber
            star_input = deserialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat")
            msk = deserialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_msk_fluxing_lst_"*lpad(adjfibindx,3,"0")*".jdat")
            subiter = star_input[.!msk]
            new_vec = map(i->map(x->x[i],subiter),1:length(subiter[1]))
            nstar = length(new_vec[1])
            new_vec[1]=1:nstar
            serialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_msked_"*lpad(adjfibindx,3,"0")*".jdat",collect(Iterators.zip(new_vec...)))
        end
    end
end

@everywhere begin
    function strip_plates(fname;nfiber=300)
        sname = split(fname,"/")
        tele = sname[end-4]
        field = sname[end-3]
        plate = sname[end-2]
        mjd = parse(Int,sname[end-1])
        sname2 = split(sname[end],"-")
        chip = sname2[2]

        f = FITS(fname)
        plugs = read(f[12],"OBJTYPE")
        close(f)
        return Iterators.zip(repeat([release_dir],nfiber),repeat([redux_ver],nfiber),repeat([tele],nfiber),repeat([field],nfiber),repeat([plate],nfiber),repeat([mjd],nfiber),1:nfiber,plugs,repeat([chip],nfiber))
    end
end
@passobj 1 workers() release_dir
@passobj 1 workers() redux_ver

if check_plates
    println("##### Using plate files to generate sky/telluric run lists #####")
    plate_flst = []
    for tele in tele_list
        if tele[1:6] =="apo25m"
            push!(plate_flst,sort(glob("*/*/*/apPlate-*",getUtahBase(release_dir,redux_ver)*"visit/$(tele)/")))
        elseif tele[1:6]=="lco25m"
            push!(plate_flst,sort(glob("*/*/*/asPlate-*",getUtahBase(release_dir,redux_ver)*"visit/$(tele)/")))
        else
            error("What telescope ($(tele)) are you talking about?")
        end 
    end
    plate_flst_all = vcat(plate_flst...)
    println("Found $(length(plate_flst_all)) plate files")
    writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_all_plates.txt",plate_flst_all,',')

    pout = @showprogress pmap(strip_plates,plate_flst_all);

    df_plate = DataFrame(
        release_dir = String[],
        redux_ver = String[],
        tele_id = String[],
        field_id = String[],
        plate_id = String[],
        mjd_id = Int[],
        fiberindx = Int[],
        fibertype = String[],
        chip_id = String[],
    );

    @showprogress for i=1:length(pout)
        for ele in pout[i]
            push!(df_plate,ele)
        end
    end

    serialize(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_df_plates.jdat",df_plate)
    df_plate_unique = unique(df_plate, [:release_dir, :redux_ver, :tele_id, :field_id, :plate_id, :mjd_id, :fiberindx, :fibertype]);
    serialize(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_df_plates_unique.jdat",df_plate_unique)

    if write_sky
        skycounts_raw = zeros(300*length(tele_list))
        skycounts_wind = zeros(300*length(tele_list))
        skycounts_wind_msked = zeros(300*length(tele_list))
        wind = 5
        for (teleind_sp,tele) in enumerate(tele_list)
            teleind = (tele[1:6] == "lco25m") ? 2 : 1
            @showprogress for fiberindx = 1:300
                adjfibindx = (teleind-1)*300 + fiberindx
                sfibindx = (length(tele_list)==2) ? adjfibindx : fiberindx
                msk = (df_plate_unique.fibertype .== "SKY") 
                msk .&= (df_plate_unique.tele_id .== tele)
                msk .&= (abs.(df_plate_unique.fiberindx .- fiberindx).==0) 
                skycounts_raw[sfibindx] = count(msk)
                msk = (df_plate_unique.fibertype .== "SKY") 
                msk .&= (df_plate_unique.tele_id .== tele)
                msk .&= (abs.(df_plate_unique.fiberindx .- fiberindx).<=wind) 
                skycounts_wind[sfibindx] = count(msk)
                skycounts_wind_msked[sfibindx] = count(msk)
                if skycounts_wind[sfibindx] > 0
                    df_sub = df_plate_unique[msk,:];
                    temp = Iterators.zip(
                        rownumber.(eachrow(df_sub)),
                        df_sub.release_dir,
                        df_sub.redux_ver,
                        df_sub.tele_id,
                        df_sub.field_id,
                        df_sub.plate_id,
                        df_sub.mjd_id,
                        df_sub.fiberindx
                    )
                    sky_input = collect(temp)
                    tout = pmap(test_apFlux,sky_input)
                    msk = (tout.!=0)
                    skycounts_wind_msked[sfibindx] = count(.!msk)

                    subiter = sky_input[.!msk]
                    new_vec = map(i->map(x->x[i],subiter),1:length(subiter[1]))
                    nstar = length(new_vec[1])
                    new_vec[1]=1:nstar
                    if skycounts_wind_msked[sfibindx]>0
                        serialize(outdir*"sky/"*"$(release_dir_n)_$(redux_ver_n)_sky_input_lst_plate_msked_"*lpad(adjfibindx,3,"0")*".jdat",collect(Iterators.zip(new_vec...)))
                    end
                end
            end
        end
        
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_skycounts_raw.txt",skycounts_raw,',')
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_skycounts_wind.txt",skycounts_wind,',')
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_skycounts_wind_msked.txt",skycounts_wind_msked,',')

        println("Sky visits per fiber: min $(minimum(skycounts_raw)), max $(maximum(skycounts_raw)), median $(median(skycounts_raw))")
        println("Sky visits per fiber with $(wind)-fiber window: min $(minimum(skycounts_wind)), max $(maximum(skycounts_wind)), median $(median(skycounts_wind))")
        println("Sky visits per fiber with $(wind)-fiber window lost to missing flux files: $(sum(skycounts_wind)-sum(skycounts_wind_msked)) , $(100*(sum(skycounts_wind)-sum(skycounts_wind_msked))/sum(skycounts_wind))%")
    end

    if write_star_plate
        ## Star from plate list (for comparison at the moment)
        starcounts_raw = zeros(300*length(tele_list))
        starcounts_msked = zeros(300*length(tele_list))
        for (teleind_sp,tele) in enumerate(tele_list)
            teleind = (tele[1:6] == "lco25m") ? 2 : 1
            @showprogress for fiberindx = 1:300
                adjfibindx = (teleind-1)*300 + fiberindx
                sfibindx = (length(tele_list)==2) ? adjfibindx : fiberindx
                msk = (df_plate_unique.fibertype .== "STAR") 
                msk .&= (df_plate_unique.tele_id .== tele)
                msk .&= (abs.(df_plate_unique.fiberindx .- fiberindx).==0) 
                starcounts_raw[sfibindx] = count(msk)
                if starcounts_raw[sfibindx] > 0
                    df_sub = df_plate_unique[msk,:];
                    temp = Iterators.zip(
                        rownumber.(eachrow(df_sub)),
                        df_sub.release_dir,
                        df_sub.redux_ver,
                        df_sub.tele_id,
                        df_sub.field_id,
                        df_sub.plate_id,
                        df_sub.mjd_id,
                        df_sub.fiberindx
                    )
                    star_input = collect(temp)
                    tout = pmap(test_apFlux,star_input)
                    msk = (tout.!=0)
                    starcounts_msked[sfibindx] = count(.!msk)

                    subiter = star_input[.!msk]
                    new_vec = map(i->map(x->x[i],subiter),1:length(subiter[1]))
                    nstar = length(new_vec[1])
                    new_vec[1]=1:nstar
                    if starcounts_msked[sfibindx]>0
                        serialize(outdir*"sky/"*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_plate_msked_"*lpad(adjfibindx,3,"0")*".jdat",collect(Iterators.zip(new_vec...)))
                    end
                end
            end
        end
        
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_starcounts_plate_raw.txt",starcounts_raw,',')
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_starcounts_plate_msked.txt",starcounts_msked,',')

        println("Star visits (plate) per fiber: min $(minimum(starcounts_raw)), max $(maximum(starcounts_raw)), median $(median(starcounts_raw))")
        println("Star visits (plate) per fiber lost to missing flux files: $(sum(starcounts_raw)-sum(starcounts_msked)) , $(100*(sum(starcounts_raw)-sum(starcounts_msked))/sum(starcounts_raw))%")
        
        run_lsts_loc = []
        run_lens_loc = []
        for telematch in tele_list
            for fiber in 1:300
                teleind = (telematch == "lco25m") ? 2 : 1
                adjfibindx = (teleind-1)*300 + fiber
                run_per_fiber = deserialize(outdir*"star/"*"$(release_dir_n)_$(redux_ver_n)_star_input_lst_msked_"*lpad(adjfibindx,3,"0")*".jdat")
                push!(run_lsts_loc,run_per_fiber)
                push!(run_lens_loc,length(run_per_fiber))
            end
        end
        run_lst_loc = vcat(run_lsts_loc...)
        
        println("Stars in visit summary not in plate files: $(length(run_lst_loc)-sum(starcounts_msked)) , $(100*(length(run_lst_loc)-sum(starcounts_msked))/length(run_lst_loc))%")
    end

    # Telluric standard star observations
    # I don't love the current input form to the Tfun script, but generate this way for now.
    # Can modify here later.
    if write_tell
        tellcounts_raw = zeros(300*length(tele_list))
        tellcounts_msked = zeros(300*length(tele_list))
        for (teleind_sp,tele) in enumerate(tele_list)
            teleind = (tele[1:6] == "lco25m") ? 2 : 1
            @showprogress for fiberindx = 1:300
                adjfibindx = (teleind-1)*300 + fiberindx
                sfibindx = (length(tele_list)==2) ? adjfibindx : fiberindx
                msk = (df_plate_unique.fibertype .== "HOT_STD") 
                msk .&= (df_plate_unique.tele_id .== tele)
                msk .&= (abs.(df_plate_unique.fiberindx .- fiberindx).==0) 
                tellcounts_raw[sfibindx] = count(msk)
                
                if tellcounts_raw[sfibindx] > 0
                    df_sub = df_plate_unique[msk,:];
                    temp = Iterators.zip(
                        rownumber.(eachrow(df_sub)),
                        df_sub.release_dir,
                        df_sub.redux_ver,
                        df_sub.tele_id,
                        df_sub.field_id,
                        df_sub.plate_id,
                        df_sub.mjd_id,
                        df_sub.fiberindx
                    )
                    tell_input = collect(temp)
                    tout = pmap(test_apFlux,tell_input)
                    msk = (tout.!=0)
                    tellcounts_msked[sfibindx] = count(.!msk)

                    subiter = tell_input[.!msk]
                    new_vec = map(i->map(x->x[i],subiter),1:length(subiter[1]))
                    nstar = length(new_vec[1])
                    new_vec[1]=1:nstar
                    if tellcounts_msked[sfibindx]>0
                        serialize(outdir*"tell/"*"$(release_dir_n)_$(redux_ver_n)_tell_input_lst_plate_msked_"*lpad(adjfibindx,3,"0")*".jdat",collect(Iterators.zip(new_vec...)))
                    end
                end
            end
        end
        
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_tellcounts_raw.txt",tellcounts_raw,',')
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_tellcounts_msked.txt",tellcounts_msked,',')

        good_tell_list = findall(tellcounts_msked.>0);
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_good_tell_list.txt",good_tell_list,',')

        println("Telluric visits (plate) per fiber: min $(minimum(tellcounts_raw)), max $(maximum(tellcounts_raw)), median $(median(tellcounts_raw))")
        println("Telluric visits (plate) per fiber lost to missing flux files: $(sum(tellcounts_raw)-sum(tellcounts_msked)) , $(100*(sum(tellcounts_raw)-sum(tellcounts_msked))/sum(tellcounts_raw))%")

        samp_lst = []
        tellcounts_wind = zeros(300*length(tele_list))
        for (teleind_sp,tele) in enumerate(tele_list)
            teleind = (tele[1:6] == "lco25m") ? 2 : 1
            for fiberindx = 1:300
                adjfibindx = (teleind-1)*300 + fiberindx
                sfibindx = (length(tele_list)==2) ? adjfibindx : fiberindx
                for wid=0:5
                    test_inds = if (teleind == 2)
                        unique(clamp.(adjfibindx .+ (-wid:wid),301,600))
                    else
                        unique(clamp.(adjfibindx .+ (-wid:wid),1,300))
                    end
                    samp_len = 0
                    for test_ind in test_inds
                        fname = outdir*"tell/"*"$(release_dir_n)_$(redux_ver_n)_tell_input_lst_plate_msked_"*lpad(test_ind,3,"0")*".jdat"
                        if isfile(fname)
                            samp_len += length(deserialize(fname))
                        end
                    end
                    if samp_len>20
                        push!(samp_lst,test_inds)
                        tellcounts_wind[sfibindx] = samp_len
                        break
                    end
                end
            end
        end
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_tellcounts_wind.txt",tellcounts_wind,',')
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_tellSamples2read.txt",samp_lst,',')

        println("Telluric visits (plate) per fiber with variable width window: min $(minimum(tellcounts_wind)), max $(maximum(tellcounts_wind)), median $(median(tellcounts_wind))")
    end
end

@everywhere begin
    function getWaveCalsFromNightRunner(runtup; cont_wid = 6)
        release_dir,redux_ver,tele,mjd = runtup
        fname = build_expPath(release_dir,redux_ver,tele,mjd)
        f = FITS(fname)
        EXPNUM = read(f[2],"NUM")
        IMAGETYP = read(f[2],"IMAGETYP")
        QRTZ = read(f[2],"QRTZ")
        THAR = read(f[2],"THAR")
        UNE = read(f[2],"UNE");
        DITHPIX = read(f[2],"DITHPIX");
        close(f)
    
        tharmsk = (IMAGETYP.=="ArcLamp") .& (THAR.==1)
        unemsk = (IMAGETYP.=="ArcLamp") .& (UNE.==1)
        fpimsk = (IMAGETYP.=="ArcLamp") .& (THAR.==0) .& (UNE.==0) .& (QRTZ.==0)
        keep_msk = tharmsk .| unemsk .| fpimsk;
        jump_ind = cumsum([0, (diff(EXPNUM[keep_msk]).>cont_wid)...]);
        
        backing_files_good = zeros(Bool,count(keep_msk))
    
        check_fpi = Iterators.zip(repeat([release_dir],count(fpimsk)),
            repeat([redux_ver],count(fpimsk)),
            repeat([tele],count(fpimsk)),
            repeat([mjd],count(fpimsk)),
            EXPNUM[fpimsk],
            )
        fpi_good = check_FPI_1d_good.(check_fpi)
        backing_files_good[fpimsk[keep_msk]].= fpi_good
    
        arc_msk = tharmsk .| unemsk
        check_arc = Iterators.zip(repeat([release_dir],count(arc_msk)),
            repeat([redux_ver],count(arc_msk)),
            repeat([tele],count(arc_msk)),
            repeat([mjd],count(arc_msk)),
            EXPNUM[arc_msk],
            )
        arc_good = check_arc_Lines_good.(check_arc)
        backing_files_good[arc_msk[keep_msk]].= arc_good
        bad_files = vcat(collect(check_fpi)[.!fpi_good],collect(check_arc)[.!arc_good])
        
        wave_run_lst = []    
        if count(keep_msk)>0
            df = DataFrame([EXPNUM[keep_msk], IMAGETYP[keep_msk], fpimsk[keep_msk], tharmsk[keep_msk], unemsk[keep_msk], DITHPIX[keep_msk], backing_files_good, jump_ind],["EXPNUM", "IMAGETYP", "FPI", "THAR", "UNE", "DITHPIX", "FILES_GOOD","CONTINDX"])
            grouped_df = groupby(df, [:DITHPIX, :CONTINDX])
    
            for group in grouped_df
                push!(wave_run_lst,wavedf2list(group,release_dir,redux_ver,tele,mjd))
            end
        end
        return wave_run_lst, bad_files, count(IMAGETYP.=="Object")
    end
    
    function wavedf2list(df_subgroup,release_dir,redux_ver,tele,mjd)
        dfs = DataFrame(df_subgroup)
        wavecallst = []
        glist = []
        while (count(dfs.THAR)>0) & (count(dfs.UNE)>0)
            sgrp = []
            # Find ThAr
            ind = findfirst(dfs.THAR)
            push!(sgrp,(release_dir,redux_ver,tele,mjd,dfs.EXPNUM[ind],dfs.DITHPIX[ind],dfs.FILES_GOOD[ind]))
            delete!(dfs, ind)
            # Find UNe
            ind = findfirst(dfs.UNE)
            push!(sgrp,(release_dir,redux_ver,tele,mjd,dfs.EXPNUM[ind],dfs.DITHPIX[ind],dfs.FILES_GOOD[ind]))
            delete!(dfs, ind)
            push!(glist,sgrp)
        end
        push!(wavecallst,glist)
        fpi_lst = []
        if count(dfs.FPI)>0
            fpi_inds = findall(dfs.FPI)
            for fpi_indx in fpi_inds
                push!(fpi_lst,(release_dir,redux_ver,tele,mjd,dfs.EXPNUM[fpi_indx],dfs.DITHPIX[fpi_indx],dfs.FILES_GOOD[fpi_indx]))
            end
        end
        push!(wavecallst,fpi_lst)
        return wavecallst
    end

    function sub_select(wave_tup)
        if length(wave_tup)>0
            fpi_grps = findall(map(x->length(x[2])!=0,wave_tup))
            for grp_indx in fpi_grps # group forces same dither
                wave_cal_lst = []
                # grab the first good arc pair
                for arc_indx in 1:length(wave_tup[grp_indx][1])
                    if (wave_tup[grp_indx][1][arc_indx][1][end] & wave_tup[grp_indx][1][arc_indx][2][end])
                        push!(wave_cal_lst,[[wave_tup[grp_indx][1][arc_indx][1][1:end-2],wave_tup[grp_indx][1][arc_indx][2][1:end-2]]])
                        break
                    end
                end
                # grab the first good fpi
                for fpi_indx in 1:length(wave_tup[grp_indx][2])
                    if wave_tup[grp_indx][2][fpi_indx][end]
                        push!(wave_cal_lst,wave_tup[grp_indx][2][fpi_indx][1:end-2])
                        break
                    end
                end
                # check if both parts of run tuple are right length, if yes return
                if (length(wave_cal_lst) == 2) && (length(wave_cal_lst[1][1]) == 2)
                    return wave_cal_lst
                end
            end
        end
        if length(wave_tup)>0
            # if no good fpi groups, loop over all groups, and only make the arc part
            for grp_indx in 1:length(wave_tup) # group forces same dither
                wave_cal_lst = []
                # grab the first good arc pair
                for arc_indx in 1:length(wave_tup[grp_indx][1])
                    if (wave_tup[grp_indx][1][arc_indx][1][end] & wave_tup[grp_indx][1][arc_indx][2][end])
                        push!(wave_cal_lst,[[wave_tup[grp_indx][1][arc_indx][1][1:end-2],wave_tup[grp_indx][1][arc_indx][2][1:end-2]]])
                        push!(wave_cal_lst,())
                        break
                    end
                end
                # check if arcpart parts of run tuple are right length, if yes return
                if (length(wave_cal_lst)>0) && (length(wave_cal_lst[1][1]) == 2)
                    return wave_cal_lst
                end
            end
        end
        wave_cal_lst = []
        return wave_cal_lst # failure, length 0
    end
    
    function check_FPI_1d_good(fpi_tup;avg_flux_cut=1e3,rad_grow=1) # could change to read file rather than just exists
        good_flag = true
        for chip in ["a","b","c"]
            good_flag &= isfile(build_framepath(fpi_tup...,chip))
            if good_flag
                ap1D_path = build_framepath(fpi_tup...,chip)
                f = FITS(ap1D_path)
                ref_img = read(f[2]);
                flg_img = read(f[4]);
                close(f)
    
                msk_bad = (flg_img .& (2^0 | 2^1 | 2^2 | 2^3 | 2^4 | 2^5 | 2^6 | 2^7 | 2^12 | 2^13 | 2^14) .!=0) #https://www.sdss4.org/dr17/irspec/apogee-bitmasks/
                msk_bad .|= grow_msk2d(flg_img .& (2^14) .!=0,rad=rad_grow)
                nan_img = copy(ref_img)
                nan_img[msk_bad].=NaN;
                good_flag &= (nansum(nan_img)./count(.!isnan.(msk_bad)) > avg_flux_cut)
            end
        end
        return good_flag
    end
    
    function check_arc_Lines_good(arc_tup) # could change to read file rather than just exists
        good_flag = true
        # for chip in ["a","b","c"]
            good_flag &= isfile(build_apLinesPath(arc_tup...))
        # end
        return good_flag
    end

    function checkIfFPIoption(wave_tup)
        fpi_option = false
        for subgrpindx in 1:length(wave_tup)
            fpi_option |= length(wave_tup[subgrpindx][2]).!=0
        end
        return fpi_option
    end

end

if check_wavecal
    println("##### Checking the files backing wavecal #####")
    wave_lst = []
    for tele in tele_list
        if tele[1:6] =="apo25m"
            tele_alias = "apogee-n"
            wave_flst = sort(glob("*/*exp.fits",getUtahBase(release_dir,redux_ver)*"exposures/$(tele_alias)/"))
            mjdlst = map(x->split(x,"/")[end-1],wave_flst);
            runlst = Iterators.zip(repeat([release_dir],length(mjdlst)),repeat([redux_ver],length(mjdlst)),repeat([tele],length(mjdlst)),mjdlst);
            push!(wave_lst,collect(runlst))
        elseif tele[1:6]=="lco25m"
            tele_alias = "apogee-s"
            wave_flst = sort(glob("*/*exp.fits",getUtahBase(release_dir,redux_ver)*"exposures/$(tele_alias)/"))
            mjdlst = map(x->split(x,"/")[end-1],wave_flst);
            runlst = Iterators.zip(repeat([release_dir],length(mjdlst)),repeat([redux_ver],length(mjdlst)),repeat([tele],length(mjdlst)),mjdlst);
            push!(wave_lst,collect(runlst))
        else
            error("What telescope ($(tele)) are you talking about?")
        end 
    end
    wave_lst_all = vcat(wave_lst...)

    pout_all = @showprogress pmap(getWaveCalsFromNightRunner,wave_lst_all);
    pout = map(x->x[1],pout_all)
    bad_files = vcat(map(x->x[2],pout_all)...)
    obj_cnts = map(x->x[3],pout_all)

    FPIpossiblevec = checkIfFPIoption.(pout)
    numFPIpossible = count(FPIpossiblevec)
    
    # subselect for current runs
    # some day, I want to just run them all (if we replace the DRP's upstream ap1D production)
    fout = sub_select.(pout);
    msklen = (length.(fout) .!=0)
    fpi_found_vec = map(x->(length(x)!=0) && (length(x[end]).!=0),fout)
    fpi_found = findall(fpi_found_vec)

    serialize(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_wavecal_input_lst.jdat",fout[msklen]);
    writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_wavecals.txt",bad_files,',')

    println("Number of missing wavecal backing files: $(length(bad_files))")
    println("MJDs without any wavecal exposures: $(count(.!msklen)) out of $(length(msklen)), $(100*count(.!msklen)/length(msklen))%")
    println("Total of $numFPIpossible nights with FPI, only $(length(fpi_found)) files found. Lost FPI for $(numFPIpossible - length(fpi_found)) nights.")
    println("Number of Object exposures that could have FPI, but don't $(sum(obj_cnts[.!fpi_found_vec .& FPIpossiblevec]))")
    println("Number of Object exposures that could have wavecal, but don't: $(sum(obj_cnts[.!msklen]))")

    mjd_vec = map(x->x[4],wave_lst_all)
    mjds_wo_cal = map(x->parse(Int,x[4]),wave_lst_all[.!msklen])
    mjds_wo_fpi = map(x->parse(Int,x[1][1][1][4]),fout[.!fpi_found_vec .& FPIpossiblevec .& msklen])
    mjds_w_science = map(x->x[7],run_lst)
    writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_mjds_w_science.txt",mjds_w_science,',')
    serialize(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_wave_obj_cnts.jdat",[mjd_vec, obj_cnts])
    mjd_wo_cal_sci_msk = mjds_wo_cal .∈ [mjds_w_science]
    mjd_wo_fpi_sci_msk = mjds_wo_fpi .∈ [mjds_w_science]
    println("Number of science MJDs without good wavecal files: $(count(mjd_wo_cal_sci_msk))")
    if count(mjd_wo_cal_sci_msk)>0
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_mjd_nowavecal.txt",mjds_wo_cal[mjd_wo_cal_sci_msk],',')
    end
    println("Number of science MJDs without good FPI files: $(count(mjd_wo_fpi_sci_msk))")
    if count(mjd_wo_fpi_sci_msk)>0
        writedlm(outdir*"summary/"*"$(release_dir_n)_$(redux_ver_n)_bad_mjd_nofpi.txt",mjds_wo_fpi[mjd_wo_fpi_sci_msk],',')
    end
end

rmprocs(workers())

# julia +1.8.2 verify_drp.jl "sdsswork/mwm" "1.2" | tee -a ../../outlists/sdsswork_1p2.log
# julia +1.8.2 verify_drp.jl "dr17" "dr17" | tee -a ../../outlists/dr17_dr17.log