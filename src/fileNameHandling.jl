### File Name Handlers

## Internal Caching Structure
function cache_skyname(tele,field,plate,mjd; cache_dir="../local_cache")
    join([cache_dir,mjd,join(["sky",tele,field,plate,mjd],"_")],"/")*".jdat"
end
#do not remake for each fiber

# choose not to cache the skyspec in injection_cache_dir (change injection tests to point cache_dir to injection_cache_dir for skyspec)
function cache_skynameSpec(tele,field,plate,mjd,fiberindx; telluric_div=false, cache_dir="../local_cache")#,inject_cache_dir="../inject_local_cache")
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    adjfibindx = (teleind-1)*300 + fiberindx
    # basedir = if tele[end]=='i'
    #     inject_cache_dir
    # else
    #     cache_dir
    # end
    basedir = cache_dir
    if telluric_div
        return join([basedir,mjd,lpad(adjfibindx,3,"0"),join(["skySpec","tellDiv",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    else
        return join([basedir,mjd,lpad(adjfibindx,3,"0"),join(["skySpec",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    end
end

function cache_starname(tele,field,plate,mjd,fiberindx; cache_dir="../local_cache",inject_cache_dir="../inject_local_cache")
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    adjfibindx = (teleind-1)*300 + fiberindx
    if tele[end]=='i'
        return join([inject_cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["star",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    else
        return join([cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["star",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    end
end

function cache_wavename(tele,mjd; cache_dir="../local_cache")
    return join([cache_dir,tele,mjd,join(["wavecal",tele,mjd],"_")],"/")*".fits"
end

function cache_waveLinesFPIname(tele,mjd,expnums,chip; cache_dir="../local_cache")
    expnum = lpad(expnums,8,"0")
    return join([cache_dir,tele,mjd,join(["waveLinesFPI",tele,mjd,expnum,chip],"_")],"/")*".fits"
end

function cache_fluxname(tele,field,plate,mjd; cache_dir="../local_cache")
    join([cache_dir,mjd,join(["fluxing",tele,field,plate,mjd],"_")],"/")*".fits"
end

## Utah SDSS Data Structure

function getUtahBase(release_dir, redux_ver)
    if occursin("dr", release_dir)
        dr_number = parse(Int, match(r"dr(\d+)", release_dir).captures[1])
        if (10 <= dr_number <= 17)
            return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/apogee/spectro/redux/$(redux_ver)/"
        elseif (18 <= dr_number)
            return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/spectro/apogee/redux/$(redux_ver)/"
        end
    end
    if release_dir[1:3] == "ipl"
        return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/spectro/apogee/redux/$(redux_ver)/"
    end
    # the last case catches the dev versions under daily/trial versions
    return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/apogee/spectro/redux/$(redux_ver)/"
end

function summary_file_by_dr(release_dir,redux_ver,dr_number,allType) # allType is either "allVisit" or "allStar"
    if dr_number==17
        return replace(getUtahBase(release_dir,redux_ver),"redux"=>"aspcap")*"synspec/$(allType)-$(redux_ver)-synspec.fits"
    elseif dr_number==16
        return replace(getUtahBase(release_dir,redux_ver),"redux"=>"aspcap")*"l33//$(allType)-$(redux_ver)-l33.fits"
    elseif dr_number==15
        return getUtahBase(release_dir,redux_ver)*"/$(allType)-l31c.2.fits"
    else
        error("No support for $(dr_number) yet, due to inhomogenous summary files")
    end
end

function build_platepath(release_dir,redux_ver,tele,field,plate,mjd,chip)
    base = getUtahBase(release_dir,redux_ver)*"visit"
    prefix = if (tele[1:6] =="apo25m")
        "apPlate"
    else
        "asPlate"
    end
    fname = "$(prefix)-$(chip)-$(plate)-$(mjd).fits"
    return join([base,tele[1:6],field,string(parse(Int,plate)),mjd,fname],"/")
end

function getFramesFromPlate(x)
    f = FITS(x)
    framenum = read(f[14],"FRAMENUM")
    close(f)
    return framenum
end

function build_framepath(release_dir,redux_ver,tele,mjd,expnums,chip)
    expnum = lpad(expnums,8,"0")
    base = getUtahBase(release_dir,redux_ver)*"exposures"
    if tele[1:6] =="apo25m"
        supfold = "apogee-n"
        fname = "ap1D-$chip-$expnum"*".fits"
        return join([base,supfold,mjd,fname],"/")
    else
        supfold = "apogee-s"
        fname = "as1D-$chip-$expnum"*".fits"
        return join([base,supfold,mjd,fname],"/")
    end
end

function build_visitpath(release_dir,redux_ver,tele,field,plate,mjd,fiberindx)
    base = getUtahBase(release_dir,redux_ver)*"visit"
    fibnum = lpad(301-fiberindx,3,"0")
    prefix = if (tele[1:6] =="apo25m")
        "apVisit"
    else
        "asVisit"
    end
    fname = "$(prefix)-$(redux_ver)-$(plate)-$(mjd)-$(fibnum).fits"
    return join([base,tele[1:6],field,string(parse(Int,plate)),mjd,fname],"/")
end

function visit2cframe(fname,tele,expnums,chip)
    expnum = lpad(expnums,8,"0")
    prefix = if (tele[1:6] =="apo25m")
        "apCframe"
    else
        "asCframe"
    end
    file = "$prefix-$chip-$expnum"*".fits"
    sname = split(fname,"/")
    sname[end] = file
    return join(sname,"/")
end

function build_expPath(release_dir,redux_ver,tele,mjd)
    base = getUtahBase(release_dir,redux_ver)*"exposures"
    fname = "$mjd"*"exp.fits"
    if tele[1:6] =="apo25m"
        supfold = "apogee-n"
        return join([base,supfold,mjd,fname],"/")
    else
        supfold = "apogee-s"
        return join([base,supfold,mjd,fname],"/")
    end
end

function build_apFPILinesPath(release_dir,redux_ver,tele,mjd,chip,expnums)
    base = getUtahBase(release_dir,redux_ver)*"cal"
    expnum = lpad(expnums,8,"0")
    if tele[1:6] =="apo25m"
        supfold = "apogee-n"
        fname = "apWaveFPI-$chip-$mjd-$expnum"*".fits"
        return join([base,supfold,"wave",fname],"/")
    else
        supfold = "apogee-s"
        fname = "asWaveFPI-$chip-$mjd-$expnum"*".fits"
        return join([base,supfold,"wave",fname],"/")
    end
end

function build_apLinesPath(release_dir,redux_ver,tele,mjd,expnums)
    base = getUtahBase(release_dir,redux_ver)*"cal"
    expnum = lpad(expnums,8,"0")
    if tele[1:6] =="apo25m"
        supfold = "apogee-n"
        fname = "apLines-$expnum"*".fits"
        return join([base,supfold,"wave",fname],"/")
    else
        supfold = "apogee-s"
        fname = "asLines-$expnum"*".fits"
        return join([base,supfold,"wave",fname],"/")
    end
end

function build_apFluxPath(release_dir,redux_ver,tele,mjd,chip,expnums)
    base = getUtahBase(release_dir,redux_ver)*"cal"
    expnum = lpad(expnums,8,"0")
    dr_number = if occursin("dr", release_dir)
        parse(Int, match(r"dr(\d+)", release_dir).captures[1])
    else
        -1
    end
    if tele[1:6] =="apo25m"
        supfold = "apogee-n"
        fname = "apFlux-$chip-$expnum"*".fits"
        if (10 <= dr_number <= 17) # are we kidding... for DR17 we did not have the supfold in there?????
            return join([base,"flux",fname],"/")
        else
            return join([base,supfold,"flux",fname],"/")
        end
    else
        supfold = "apogee-s"
        fname = "asFlux-$chip-$expnum"*".fits"
        if (10 <= dr_number <= 17)
            return join([base,"flux",fname],"/")
        else
            return join([base,supfold,"flux",fname],"/")
        end
    end
end

# a lot of replicated code with the verify_drp.jl apFlux checker... would be good to consolidate
function build_apFluxPaths(release_dir,redux_ver,tele,field,plate,mjd)
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
        error("No DomeFlats taken")
    elseif !isnothing(expIndex_before) && (CARTID[expIndex_before] == cartVisit)
        expIndex_before
    elseif !isnothing(expIndex_after) && (CARTID[expIndex_after] == cartVisit)
        expIndex_after
    else
        error("No DomeFlats taken with the right CARTID")
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
                error("No apFlux files exist for the whole (FPS era) night")
            end
        end
    end

    expectDome = EXPNUM[expIndex]

    fluxFiles = []
    for chip in ["a","b","c"]
        fluxFile = build_apFluxPath(release_dir,redux_ver,tele,mjd,chip,expectDome)
        push!(fluxFiles,fluxFile)
    end
    return fluxFiles, expectDome, cartVisit
end