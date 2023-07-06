## File Name Handlers

function cache_skyname(intup; cache_dir="../local_cache")
    (tele,field,plate,mjd,file,plateFile,fiberindx) = intup
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    adjfibindx = (teleind-1)*300 + fiberindx
    join([cache_dir,mjd,join(["sky",tele,field,plate,mjd],"_")],"/")*".jdat"
end

function cache_skynameSpec(intup; cache_dir="../local_cache",inject_cache_dir="../inject_local_cache")
    (tele,field,plate,mjd,file,plateFile,fiberindx) = intup
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    adjfibindx = (teleind-1)*300 + fiberindx
    if tele[end]=='i'
        return join([inject_cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["skySpec",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    else
        return join([cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["skySpec",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    end
end

function cache_starname(intup; cache_dir="../local_cache",inject_cache_dir="../inject_local_cache")
    (tele,field,plate,mjd,file,plateFile,fiberindx) = intup
    teleind = (tele[1:6] == "lco25m") ? 2 : 1
    adjfibindx = (teleind-1)*300 + fiberindx
    if tele[end]=='i'
        return join([inject_cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["star",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    else
        return join([cache_dir,mjd,lpad(adjfibindx,3,"0"),join(["star",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
    end
end

function getFramesFromPlate(x)
    f = FITS(x)
    framenum = read(f[14],"FRAMENUM")
    close(f)
    return framenum
end

function build_framepath(tele,mjd,imid,chip)
    imids = lpad(imid,8,"0")
    base = "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/exposures"
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

function build_visitpath(intup)
    (tele,field,plate,mjd,file) = intup
    base = "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit"
    return join([base,tele[1:6],field,string(parse(Int,plate)),mjd,file],"/")
end

function visit2cframe(fname,tele,imid,chip)
    imids = lpad(imid,8,"0")
    if tele[1:6] =="apo25m"
        file = "apCframe-$chip-$imids"*".fits"
        sname = split(fname,"/")
        sname[end] = file
        return join(sname,"/")
    else
        file = "asCframe-$chip-$imids"*".fits"
        sname = split(fname,"/")
        sname[end] = file
        return join(sname,"/")
    end
end

function plate2visit(plate_path,fiberindx)
    fibnum = lpad(301-fiberindx,3,"0")
    sname = split(plate_path,"/")
    tele, field, plate, mjd = sname[end-4], sname[end-3], sname[end-2], sname[end-1]
    file = replace(sname[end],"apPlate"=>"apVisit")
    file = replace(file,"asPlate"=>"asVisit")
    
    file = replace(file,"-a-"=>"-dr17-")
    file = replace(file,"-b-"=>"-dr17-")
    file = replace(file,"-c-"=>"-dr17-")
    
    file = replace(file,".fits"=>"-"*fibnum*".fits")
    return file
end
