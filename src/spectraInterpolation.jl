## Spectra Interpolation Module

function find_yinx(x,y)
    lx = length(x)
    ly = length(y)
    out = zeros(Int,ly)
    @inbounds for i=1:length(y)
        ind = searchsortedfirst(x,y[i])
        if ind == 1
            out[i] = ind
        elseif ind > lx
            out[i] = ind-1
        else
            if abs(x[ind]-y[i]) < abs(x[ind-1]-y[i])
                out[i] = ind
            else
                out[i] = ind-1
            end
        end
    end
    return out
end

function returnWeights(modCoordAll,targVal,cindx;kernsize=4,kerntrim=true)
    modlen = length(modCoordAll)
    pscale = minimum(abs.(diff(modCoordAll[maximum([1,(cindx-1)]):minimum([(cindx+1),modlen])])))
    offset = (modCoordAll[cindx].-targVal)/pscale

    indvec = (-kernsize:kernsize) .+ cindx
    offvec = (-kernsize:kernsize) .+ offset
    msk = (1 .<= indvec .<= modlen) # within bounds range
    msk .&= (-kernsize .<= offvec .<= kernsize) # within kernel bounds

    if kerntrim
        if count(msk) >= 2*kernsize
            indvecr = indvec[msk]
            wvec = Interpolations.lanczos.(offvec[msk],kernsize)
            return indvecr, wvec
        else
            return zeros(Int,2*kernsize), NaN*ones(2*kernsize)
        end
    else
        indvecr = indvec[msk]
        wvec = Interpolations.lanczos.(offvec[msk],kernsize)
        return indvecr, wvec
    end
end

function returnWeights_inv(obsCoordall,obsBitMsk,targVal,cindx;kernsize=4)
    obslen = length(obsCoordall)
    pscale = minimum(abs.(diff(obsCoordall[maximum([1,(cindx-1)]):minimum([(cindx+1),obslen])])))
    offset = (obsCoordall[cindx].-targVal)/pscale
    cbit = obsBitMsk[cindx]
    if (cbit .& 2^1)!=0
        cchip = 1
    elseif (cbit .& 2^2)!=0
        cchip = 2
    elseif (cbit .& 2^3)!=0
        cchip = 3
    else
        cchip = 0
        # println("NO CHIP?")
    end

    indvec = (-kernsize:kernsize) .+ cindx
    offvec = (-kernsize:kernsize) .+ offset
    msk = (1 .<= indvec .<= obslen) # within bounds range
    msk .&= (-kernsize .< offvec .< kernsize) # within kernel bounds
    indvecr = indvec[msk]

    mskb = ((obsBitMsk[indvecr] .& 2^cchip).!=0) # same chip mask
    mskb .&= ((obsBitMsk[indvecr] .& 2^4).==0) #bad pix mask
    indvecrr = indvecr[mskb]

    wvec = Interpolations.lanczos.(offvec[msk][mskb],kernsize)
    return indvecrr, wvec
end

function generateInterpMatrix_sparse_inv(waveobs,obsBitMsk,wavemod)
    obslen = length(waveobs)
    modlen = length(wavemod)
    cindx = find_yinx(waveobs,wavemod)
    row, col, val = Int[], Int[], Float64[]
    for (modind, modval) in enumerate(wavemod)
        indxvec, wvec = returnWeights_inv(waveobs,obsBitMsk,modval,cindx[modind])
        wvec ./= sum(wvec)
        if length(indxvec) > 7
            push!(row, (modind.*ones(Int,length(indxvec)))...)
            push!(col, indxvec...)
            push!(val, wvec...)
        end
    end
    return sparse(row,col,val,modlen,obslen)
end

function generateInterpMatrix_sparse(waveobs,wavemod;kernsize=4,kerntrim=true)
    obslen = length(waveobs)
    modlen = length(wavemod)
    cindx = find_yinx(wavemod,waveobs)
    row, col, val = Int[], Int[], Float64[]
    for (obsin, obsval) in enumerate(waveobs)
        indxvec, wvec = returnWeights(wavemod,obsval,cindx[obsin],kernsize=kernsize,kerntrim=kerntrim)
        if !isnan(wvec[1])
            wvec ./= sum(wvec)
            push!(row, indxvec...)
            push!(col, (obsin.*ones(Int,length(indxvec)))...)
            push!(val, wvec...)
        end
    end
    return sparse(row,col,val,modlen,obslen)
end

function generateSelectMatrix_sparse(waveobs,obsBitMsk)
    obslen = length(waveobs)
    goodpix = ((obsBitMsk .& 2^4).==0)
    goodlen = count(goodpix)
    row, col, val = Int[], Int[], Float64[]
    gindx = 1
    for ind=1:obslen
        if goodpix[ind]
            push!(row, ind)
            push!(col, gindx)
            push!(val, 1)
            gindx +=1
        end
    end
    return sparse(row,col,val,obslen,goodlen)
end
