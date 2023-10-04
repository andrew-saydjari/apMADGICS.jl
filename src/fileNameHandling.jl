### File Name Handlers

## Internal Caching Structure
function cache_skyname(tele,field,plate,mjd; cache_dir="../local_cache")
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    join([cache_dir,mjd,join(["sky",tele,field,plate,mjd],"_")],"/")*".jdat"
end
#do not remake for each fiber

function cache_skynameSpec(tele,field,plate,mjd,fiberindx; cache_dir="../local_cache")
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    adjfibindx = (teleind-1)*300 + fiberindx
    return join([cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["skySpec",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
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

## Utah SDSS Data Structure

function getUtahBase(release_dir, redux_ver)
    if occursin("dr", release_dir)
        dr_number = parse(Int, match(r"dr(\d+)", redux_ver).captures[1])
        if (10 <= dr_number <= 17)
            return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/apogee/spectro/redux/$(redux_ver)/"
        end
    end
    if release_dir[1:3] == "ipl"
        return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/spectro/apogee/redux/$(redux_ver)/"
    end
    return "/uufs/chpc.utah.edu/common/home/sdss/$(release_dir)/apogee/spectro/redux/$(redux_ver)/"
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

function build_framepath(release_dir,redux_ver,tele,mjd,imid,chip)
    imids = lpad(imid,8,"0")
    base = getUtahBase(release_dir,redux_ver)*"exposures"
    if tele[1:6] =="apo25m"
        supfold = "apogee-n"
        fname = "ap1D-$chip-$imids"*".fits"
        return join([base,supfold,mjd,fname],"/")
    else
        supfold = "apogee-s"
        fname = "as1D-$chip-$imids"*".fits"
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

function visit2cframe(fname,tele,imid,chip)
    imids = lpad(imid,8,"0")
    prefix = if (tele[1:6] =="apo25m")
        "apCframe"
    else
        "asCframe"
    end
    file = "$prefix-$chip-$imids"*".fits"
    sname = split(fname,"/")
    sname[end] = file
    return join(sname,"/")
end
