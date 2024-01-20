## Ingest Module
# minor rewrite on the horizon

using AstroTime

function getAndWrite_fluxing(release_dir,redux_ver,tele,field,plate,mjd; cache_dir="../local_cache")
    flux_paths, domeflat_expid = build_apFluxPaths(release_dir,redux_ver,tele,field,plate,mjd)
    fluxingcache = cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir)

    hdr = FITSHeader(["pipeline","git_branch","git_commit","domeflat_expid"],["apMADGICS.jl",git_branch,git_commit,string(domeflat_expid)],["","","",""])
    h = FITS(fluxingcache,"w")
    write(h,[0],header=hdr,name="header_only")
    for (chipind,chip) in enumerate(["a","b","c"])
        flux_path = flux_paths[chipind]
        f = FITS(flux_path)
        thrpt = read(f[3])
        close(f)
        write(h,thrpt,name=chip)
    end
    close(h)
end

function getSky4visit(release_dir,redux_ver,tele,field,plate,mjd,fiberindx; caching=false,cache_dir="../local_cache")
    
    ### Find all of the Sky Fibers
    vname = build_visitpath(release_dir,redux_ver,tele,field,plate,mjd,fiberindx)
    # do we want to move away from relying on the apPlateFile? 
    # idk, seems fairly robust
    plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
    frame_lst = getFramesFromPlate(plateFile)

    fname = visit2cframe(vname,tele,frame_lst[1],"a")
    f = FITS(fname)
    objtype = read(f[12],"OBJTYPE")
    fiberids = read(f[12],"FIBERID");
    close(f)
    skyinds = findall((objtype.=="SKY") .& (fiberids.!=(301-fiberindx)));

    ### Decompose all of the Sky Fibers
    outcont = zeros(length(wavetarg),length(skyinds));
    for (findx,fiberind) in enumerate(skyinds)
        try
            skycacheSpec = cache_skynameSpec(tele,field,plate,mjd,fiberind,cache_dir=cache_dir)
            if (isfile(skycacheSpec) & caching)
                fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skycacheSpec)
                starscale,framecnts,varoffset,varflux = metaexport
            else
                fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberind,cache_dir=cache_dir)
                starscale,framecnts,varoffset,varflux = metaexport
                if caching
                    dirName = splitdir(skycacheSpec)[1]
                    if !ispath(dirName)
                        mkpath(dirName)
                    end
                    serialize(skycacheSpec,[fvec, fvarvec, cntvec, chipmidtimes, metaexport])
                end
            end

            simplemsk = (cntvec.==maximum(cntvec)) .& skymsk;
            contvec = sky_decomp(fvec, fvarvec, simplemsk)
            # do we want to save the other components to disk? I am not sure we do.
            outcont[:,findx] .= contvec
        catch
            # we should figure out how often this happens and why
            @warn "Bad Sky: ($tele,$field,$plate,$mjd,$fiberind)"
            outcont[:,findx] .= 0
        end
    end

    msk = (dropdims(nansum(outcont,1),dims=1)) .> 0
    meanLocSky = dropdims(nanmean(outcont[:,msk],2),dims=2);

    VLocSky = (outcont[:,msk].-meanLocSky)./sqrt(count(msk));
    return meanLocSky, VLocSky
end

function sky_decomp(outvec,outvar,simplemsk)   
    ## Select data for use (might want to handle mean more generally)
    Xd_obs = outvec[simplemsk];

    ## Set up residuals prior
    A = Diagonal(outvar[simplemsk]);
    Ainv = Diagonal(1 ./outvar[simplemsk]);

    ## Set up priors
    V_skyline_c = V_skyline
    V_skyline_r = V_skyline_c[simplemsk,:]
    V_skycont_c = V_skycont
    V_skycont_r = V_skycont_c[simplemsk,:]
    
    # Compute sky line/continuum separation
    Vcomb = hcat(V_skyline_r,V_skycont_r);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, (V_skycont_r, ), (V_skycont_c, ))

    return x_comp_lst[1]
end

function stack_out(release_dir,redux_ver,tele,field,plate,mjd,fiberindx; varoffset=16.6, cache_dir="../local_cache")

    # # hardcoded flux-dep variance correction (empitical IPC + LSF correction)
    # power 2 model
    (p, c) = if (tele[1:6] == "apo25m")
        (2.0, 1.7e-2)    
    else
        (2.0, 3.4e-2)
    end

    plateFile = build_platepath(release_dir,redux_ver,tele,field,plate,mjd,"a")
    frame_lst = getFramesFromPlate(plateFile)
    
    #make a dictionary of values for chip a,b,c from fluxing fits
    thrptDict = Dict{String,Float64}()
    f = FITS(cache_fluxname(tele,field,plate,mjd; cache_dir=cache_dir))
    for chip in ["a","b","c"]
        thrpt = read(f[chip],fiberindx)
        thrptDict[chip] = thrpt
    end
    close(f)
    
    fill!(outvec,0)
    fill!(outvar,0)
    fill!(cntvec,0)
    time_lsts = [[],[],[]]
    for imid in frame_lst
        fill!(Xd_stack,0)
        fill!(Xd_std_stack,0)
        fill!(waveobs_stack,0)
        fill!(pixmsk_stack,0)
        fill!(telluric_stack,0)
        fill!(fullBit,0)
        for (chipind,chip) in enumerate(["c","b","a"]) #needs to be c,b,a for chip ind to be right
            fname = build_framepath(release_dir,redux_ver,tele,mjd,imid,chip)
            f = FITS(fname)
            hdr = read_header(f[1])
            midtime = modified_julian(TAIEpoch(hdr["DATE-OBS"]))+(hdr["EXPTIME"]/2/3600/24)days #TAI or UTC?
            push!(time_lsts[chipind],AstroTime.value(midtime))
            Xd = read(f[2],:,fiberindx)./thrptDict[chip];
            Xd_stack[(1:2048).+(chipind-1)*2048] .= Xd[end:-1:1]
            Xd_std = read(f[3],:,fiberindx)./thrptDict[chip];
            Xd_std_stack[(1:2048).+(chipind-1)*2048] .= Xd_std[end:-1:1]
            pixmsk = read(f[4],:,fiberindx);
            pixmsk_stack[(1:2048).+(chipind-1)*2048] .= pixmsk[end:-1:1]
            waveobsa = read(f[5],:,fiberindx);
            waveobs_stack[(1:2048).+(chipind-1)*2048] .= waveobsa[end:-1:1]
            fullBit[(1:2048).+(chipind-1)*2048] .+= 2^chipind
            close(f)
        end
        fullBit[((pixmsk_stack .& 2^0).!=0)] .+= 2^4 # call pixmask bit 0 bad
        fullBit[fullBit.==0] .+= 2^4 # call chip gaps bad for alt space

        goodpix = ((pixmsk_stack .& 2^0).==0) .& ((fullBit .& 2^4).==0)

        obsBit = fullBit[goodpix]
        Xd_obs = Xd_stack[goodpix]
        Xd_std_obs = Xd_std_stack[goodpix];
        waveobs = waveobs_stack[goodpix];
        pixindx = (1:length(waveobs_stack))[goodpix]

        Rinv = generateInterpMatrix_sparse_inv(waveobs,obsBit,wavetarg,pixindx)
        normvec = dropdims(sum(Rinv,dims=2),dims=2)
        msk_inter = (normvec.!=0)

        fullvec = Rinv*Xd_obs
        fullvec[.!msk_inter] .= 0
        varvec =  (Rinv.^2)*(Xd_std_obs.^2)
        varvec[.!msk_inter] .= 0

        outvec .+= fullvec
        outvar .+= varvec
        cntvec .+= msk_inter
    end
    framecnts = maximum(cntvec)
    outvec./=framecnts
    outvar./=(framecnts^2)
    
    simplemsk = (cntvec.==framecnts)
    starscale = if count(simplemsk .& (.!isnan.(outvec)))==0
        NaN
    else
        abs(nanmedian(outvec[simplemsk]))
    end
    # this is a systematic correction to the variance (~ 4ADU to the uncertainties) to prevent chi2 versus frame number trends
    outvar .+= varoffset
    # this is an empirical LSF and IPC correction to the variance
    outvar .+= (c^2*starscale^p)

    goodframeIndx = length.(time_lsts).!=0
    chipmidtimes = zeros(3)
    chipmidtimes[goodframeIndx] .= mean.(time_lsts[goodframeIndx]) #consider making this flux weighted (need to worry about skyline variance driving it)
    chipmidtimes[.!goodframeIndx] .= NaN
    metaexport = (starscale,framecnts,varoffset,(c^2*starscale^p),thrptDict["a"],thrptDict["b"],thrptDict["c"])
    return outvec, outvar, cntvec, chipmidtimes, metaexport
end
