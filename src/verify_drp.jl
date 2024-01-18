##### Verify APOGEE DRP Products from Release are Production Ready

import Pkg
using InteractiveUtils; versioninfo()
Pkg.activate("../"); Pkg.instantiate(); Pkg.precompile()

using FITSIO, StatsBase, ProgressMeter, Distributed, Serialization, Glob

src_dir = abspath("./")
include(src_dir*"/fileNameHandling.jl")

outdir = "../../../2024_01_09/outlists/"
if !ispath(outdir)
    mkpath(outdir)
end

release_dir = "sdsswork/mwm"
redux_ver = "1.2"

check_ap1d = false
check_apCframes = true

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

    using FITSIO, StatsBase, ProgressMeter, Distributed, Serialization, Glob

    src_dir = abspath("./")
    include(src_dir*"/fileNameHandling.jl")
end

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
                serialize(outdir*"star_input_lst_"*lpad(adjfibindx,3,"0")*".jdat",collect(subiter)[outcheck])
                if count(.!outcheck)>0
                    serialize(outdir*"star_bad_lst_"*lpad(adjfibindx,3,"0")*".jdat",collect(subiter)[.!outcheck])
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
        writedlm(outdir*"$(release_dir)_$(redux_ver)_bad_ap1D.txt",badfilevec,',')
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
        writedlm(outdir*"$(release_dir)_$(redux_ver)_bad_apCframe.txt",fnames[.!ocheck],',')
    end
end

rmprocs(workers())

