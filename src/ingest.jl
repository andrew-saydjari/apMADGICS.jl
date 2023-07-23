## Ingest Module
# minor rewrite on the horizon

using AstroTime

function getSky4visit(intup; caching=false, inject_cache_dir="../inject_local_cache",cache_dir="../local_cache")
    (tele,field,plate,mjd,file,plateFile,fiber) = intup
    
    # Find all of the Sky Fibers
    vname = build_visitpath((tele,field,plate,mjd,file))
    
    # do we want to move away from relying on the apPlateFile? 
    # idk, seems fairly robust
    frame_lst = getFramesFromPlate(plateFile)

    fname = visit2cframe(vname,tele,frame_lst[1],"a")
    f = FITS(fname)
    objtype = read(f[12],"OBJTYPE")
    fiberids = read(f[12],"FIBERID");
    close(f)
    skyinds = findall((objtype.=="SKY") .& (fiberids.!=(301-fiber)));

    outcont = zeros(length(wavetarg),length(skyinds));
    # Decompose all of the Sky Fibers
    for (findx,fiberind) in enumerate(skyinds)
        tintup = (tele,field,plate,mjd,file,plateFile,fiberind)
        try
            skycacheSpec = cache_skynameSpec(tintup,inject_cache_dir=inject_cache_dir,cache_dir=cache_dir)
            if (isfile(skycacheSpec) & caching)
                fvec, fvarvec, cntvec, chipmidtimes, metaexport = deserialize(skycacheSpec)
                starscale,framecnts,varoffset,varflux = metaexport
            else
                fvec, fvarvec, cntvec, chipmidtimes, metaexport = stack_out(tintup)
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

function stack_out(intup; varoffset=16.6)
    (tele,field,plate,mjd,file,plateFile,fiber) = intup

    # # hardcoded flux-dep variance correction (empitical IPC + LSF correction)
    # power 2 model
    (p, c) = if (tele[1:6] == "apo25m")
        (2.0, 1.7e-2)    
    else
        (2.0, 3.4e-2)
    end

    frame_lst = getFramesFromPlate(plateFile)

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
            fname = build_framepath(tele,mjd,imid,chip)
            f = FITS(fname)
            hdr = read_header(f[1])
            midtime = modified_julian(TAIEpoch(hdr["DATE-OBS"]))+(hdr["EXPTIME"]/2/3600/24)days #TAI or UTC?
            push!(time_lsts[chipind],AstroTime.value(midtime))
            Xd = read(f[2],:,fiber);
            Xd_stack[(1:2048).+(chipind-1)*2048] .= Xd[end:-1:1]
            Xd_std = read(f[3],:,fiber);
            Xd_std_stack[(1:2048).+(chipind-1)*2048] .= Xd_std[end:-1:1]
            pixmsk = read(f[4],:,fiber);
            pixmsk_stack[(1:2048).+(chipind-1)*2048] .= pixmsk[end:-1:1]
            waveobsa = read(f[5],:,fiber);
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
    ## Having made this change 05/28/2023, we now need to be careful to not divide later in any of the prior scripts
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
    metaexport = (starscale,framecnts,varoffset,(c^2*starscale^p))
    return outvec, outvar, cntvec, chipmidtimes, metaexport
end
