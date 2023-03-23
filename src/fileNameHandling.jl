## File Name Handlers

function cache_skyname(intup; cache_dir="../local_cache")
    (tele,field,plate,mjd,file,plateFile,fiberindx) = intup
    join([cache_dir,join(["sky",tele,field,plate,mjd],"_")],"/")*".jdat"
end

function cache_starname(intup; cache_dir="../local_cache")
    (tele,field,plate,mjd,file,plateFile,fiberindx) = intup
    join([cache_dir,join(["star",tele,field,plate,mjd,fiberindx],"_")],"/")*".jdat"
end

function getFramesFromPlate(x)
    f = FITS(x)
    framenum = read(f[14],"FRAMENUM")
    close(f)
    return framenum
end

function build_framepath(mjd,imid,chip)
    imids = lpad(imid,8,"0")
    base = "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/exposures/apogee-n"
    fname = "ap1D-$chip-$imids"*".fits"
    return join([base,mjd,fname],"/")
end

function build_visitpath(intup)
    (tele,field,plate,mjd,file) = intup
    base = "/uufs/chpc.utah.edu/common/home/sdss/dr17/apogee/spectro/redux/dr17/visit"
    return join([base,tele,field,string(parse(Int,plate)),mjd,file],"/")
end

function visit2cframe(fname,imid,chip)
    imids = lpad(imid,8,"0")
    file = "apCframe-$chip-$imids"*".fits"
    sname = split(fname,"/")
    sname[end] = file
    return join(sname,"/")
end

function platepath2intuple(plate_path)
    sname = split(plate_path,"/")
    field, plate, mjd = sname[end-3], sname[end-2], sname[end-1]
    file = replace(replace(replace(replace(replace(sname[end],"apPlate"=>"apVisit"),"-a-"=>"-dr17-"),"-b-"=>"-dr17-"),"-c-"=>"-dr17-"),".fits"=>"-006.fits")
    return ("apo25m",field,plate,mjd,file,295)
end
