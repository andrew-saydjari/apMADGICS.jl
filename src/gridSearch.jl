## GridSearch Module

function shift_range(inrng,shift)
    start, stop, step1 = inrng[1], inrng[end], step(inrng)
    return range(start+shift,stop+shift,step=step1)
end

function shift_trim_range(inrng,shift; minv=4//10, maxv=4)
    start, stop, step1 = inrng[1], inrng[end], step(inrng)
    return range(maximum((start+shift,minv)),minimum((stop+shift,maxv)),step=step1)
end

function quadratic_interp3(vrng,chi2,minInd;step=1)
    lval = length(vrng)
    subgrid = minInd .+ (-step:step)
    x1, x2, x3 = vrng[subgrid]
    c1, c2, c3 = chi2[subgrid]
    
    x1x2 = x1-x2
    x1x3 = x1-x3
    x2x3 = x2-x3
    k1 = c1/(x1x2*x1x3)
    k2 = c2/(-x1x2*x2x3)
    k3 = c3/(x2x3*x1x3)
    xopt = k1*(x2+x3)+k2*(x1+x3)+k3*(x2+x1)
    xopt /= (2*(k1+k2+k3))
    copt = k1*(xopt-x2)*(xopt-x3)+k2*(xopt-x1)*(xopt-x3)+k3*(xopt-x2)*(xopt-x1)
    if copt > c2 # this should not be possible and should never be thrown
        flag = 2^0
        xoptr = x2
        coptr = c2
    else
        flag = 0
        xoptr = xopt
        coptr = copt
    end
    return coptr, xoptr, flag
end

function quadratic_interp3_2d(vrngs,chi2,minInd;step1=1,step2=1)
    vrng1, vrng2 = vrngs
    lval1 = length(vrng1)
    lval2 = length(vrng2)
    subgrid1 = minInd[1] .+ (-step1:step1)
    subgrid2 = minInd[2] .+ (-step2:step2)

    if (lval1 == 1) & (lval2 != 1)
        copt, xopt, flag = quadratic_interp3(vrng2,view(chi2,1,:),minInd[2],step=step2)
        return copt, (vrng1[minInd[1]], xopt), flag
    elseif (lval1 != 1) & (lval2 == 1)
        copt, xopt, flag = quadratic_interp3(vrng1,view(chi2,:,1),minInd[1],step=step1)
        return copt, (xopt, vrng2[minInd[2]]), flag
    elseif (lval1 != 1) & (lval2 != 1)
        locgrid = chi2[subgrid1,subgrid2]
        x1, x2, x3 = vrng1[subgrid1]
        dx = (x2-x1)
        if !((x3-x2) ≈ dx)
            @warn "Non uniform grid error estimates not implemented"
        end
        y1, y2, y3 = vrng2[subgrid2]
        dy = (y2-y1)
        if !((y3-y2) ≈ dy)
            @warn "Non uniform grid error estimates not implemented"
        end
        Hxx = (locgrid[3,2]-2*locgrid[2,2]+locgrid[1,2])/(dx^2)
        Hyy = (locgrid[2,3]-2*locgrid[2,2]+locgrid[2,1])/(dy^2)
        Hxy = (locgrid[3,3]-locgrid[3,1]-locgrid[1,3]+locgrid[1,1])/(4*dx*dy)
        det0 = Hxx*Hyy-Hxy^2
        errx= Hyy/det0
        erry= Hxx/det0
        errxy= -Hxy/det0

        Gx = (locgrid[3,2]-locgrid[1,2])/(2*dx)
        Gy = (locgrid[2,3]-locgrid[2,1])/(2*dy)

        x0 = [x2, y2]
        G = [Gx; Gy]
        H = [Hxx Hxy; Hxy Hyy]
        Hinv = [errx errxy; errxy erry]

        minVal_interp =  x0 .- (Hinv*G)
        function f(x)
            dx = (x-x0)
            return locgrid[2,2] + G'*dx + (dx'*(H*dx))/2
        end
        minChi_interp = f(minVal_interp)

        # In 1d, this is definitely not possible
        # in 2d, this flag might be important
        if minChi_interp > locgrid[2,2]
            flag = 2^0
            minVal_opt = minInd
            copt = locgrid[2,2]
        else
            flag = 0
            minVal_opt = minVal_interp
            copt = minChi_interp
        end

        return copt, minVal_opt, flag
    else
        @warn "I can't interpolate a 0d scan..."
    end
end

function sampler_1d_dense(chi2_fun,valrng)
    lval = length(valrng)
    chi2out = zeros(lval)
    for (ind,val) in enumerate(valrng)
        chi2out[ind] = chi2_fun(val)
    end
    minChi, minInd = findmin(chi2out)
    flag = (2^1)*((minInd == 1) | (minInd == lval)) # edge of grid bit
    if flag !=0
        return ((valrng[minInd], minChi, valrng[minInd], minChi, minInd, flag), valrng, chi2out)
    else
        minChi_interp, minVal_interp, flag = quadratic_interp3(valrng,chi2out,minInd,step=1)
        return ((minVal_interp, minChi_interp, valrng[minInd], minChi, minInd, flag), valrng, chi2out)
    end
end

function sampler_1d_dense_var(chi2_fun,valrng;stepx=1)
    lval = length(valrng)
    chi2out = zeros(lval)
    for (ind,val) in enumerate(valrng)
        chi2out[ind] = chi2_fun(val)
    end
    minChi, minInd = findmin(chi2out)
    flag = (2^1)*((minInd == 1) | (minInd == lval)) # edge of grid bit
    if flag !=0
        return ((valrng[minInd], minChi, valrng[minInd], minChi, minInd, flag, NaN), valrng, chi2out)
    else
        minChi_interp, minVal_interp, flag = quadratic_interp3(valrng,chi2out,minInd,step=1)
        varx = err1d(valrng,chi2out,minInd;stepx=stepx)
        if isnan(varx)
            flag += 2^3 #var off grid
        end
        if varx<=0
            flag += 2^5 #bad error curvature
        end
        return ((minVal_interp, minChi_interp, valrng[minInd], minChi, minInd, flag, varx), valrng, chi2out)
    end
end

function sampler_1d_hierarchy_var(chi2_fun,lvltup;minres=1//10,stepx=1)
    lvllen = length(lvltup)
    out = []
    
    global_minVal_interp = 0
    global_minChi_interp = 0
    global_minVal = 0
    global_minChi = 0
    global_minInd = 1
    global_flag = 0
    for (lvlind,lvl) in enumerate(lvltup)
        if lvlind == 1
            uprng = lvl
        else
            uprng = shift_range(lvl,round_step(global_minVal_interp,minres))
        end
        ((minVal_interp, minChi_interp, minVal, minChi, minInd, flag), valrng, chi2out) = sampler_1d_dense(chi2_fun,uprng)
        push!(out,((minVal_interp, minChi_interp, minVal, minChi, minInd, flag), valrng, chi2out))
        if minChi < global_minChi # this should not be the interpolated value, it should be the min measured chi2 value
            global_minVal_interp, global_minChi_interp, global_minVal, global_minChi, global_minInd, global_flag = minVal_interp, minChi_interp, minVal, minChi, minInd, flag
        end
        if lvlind == lvllen
            varx = err1d(lvl,chi2out,minInd,stepx=stepx)
            if isnan(varx)
                global_flag += 2^3 #var off grid
            end
            if varx<=0
                global_flag += 2^5 #bad error curvature
            end
            return ((global_minVal_interp, global_minChi_interp, global_minVal, global_minChi, global_minInd, global_flag, varx), out)
        end
    end
end

function sampler_2d_dense(chi2_fun,valrngtup)
    valrng1, valrng2 = valrngtup
    lval1 = length(valrng1)
    lval2 = length(valrng2)
    chi2out = zeros(lval1,lval2)
    for (ind2,val2) in enumerate(valrng2)
        for (ind1,val1) in enumerate(valrng1)
            chi2out[ind1,ind2] = chi2_fun((val1,val2))
        end
    end
    minChi, minInd = findmin(chi2out)
    # now we have a gridx and gridy edge bit (might need to update bit identities)
    flag = (2^1)*((minInd[1] == 1) | (minInd[1] == lval1))*(lval1!=1) + (2^2)*((minInd[2] == 1) | (minInd[2] == lval2))*(lval2!=1)
    if flag !=0
        return (((valrng1[minInd[1]], valrng2[minInd[2]]), minChi, (valrng1[minInd[1]], valrng2[minInd[2]]), minChi, minInd, flag), valrngtup, chi2out)
    else
        minChi_interp, minVal_interp, flag = quadratic_interp3_2d(valrngtup,chi2out,minInd;step1=1,step2=1)
        return (((minVal_interp[1], minVal_interp[2]), minChi_interp, (valrng1[minInd[1]], valrng2[minInd[2]]), minChi, minInd, flag), valrngtup, chi2out)
    end
    # we may want 1d quadratic interpolations here on if cases
end

function sampler_2d_dense_var(chi2_fun,valrngtup;step1=1,step2=1)
    valrng1, valrng2 = valrngtup
    lval1 = length(valrng1)
    lval2 = length(valrng2)
    chi2out = zeros(lval1,lval2)
    for (ind2,val2) in enumerate(valrng2)
        for (ind1,val1) in enumerate(valrng1)
            chi2out[ind1,ind2] = chi2_fun((val1,val2))
        end
    end
    minChi, minInd = findmin(chi2out)
    # now we have a gridx and gridy edge bit (might need to update bit identities)
    flag = (2^1)*((minInd[1] == 1) | (minInd[1] == lval1)) + (2^2)*((minInd[2] == 1) | (minInd[2] == lval2))
    if flag !=0
        return (((valrng1[minInd[1]], valrng2[minInd[2]]), minChi, (valrng1[minInd[1]], valrng2[minInd[2]]), minChi, minInd, flag, NaN, NaN, NaN, NaN, NaN), valrngtup, chi2out)
    else
        minChi_interp, minVal_interp, flag = quadratic_interp3_2d(valrngtup,chi2out,minInd;step1=1,step2=1)
        # we may want 1d quadratic interpolations here on if cases

        # one could do 1d error bars on edge of grid problems with careful if cases
        varx, vary, varxy, varxx, varyy = err2d((valrng1,valrng2),chi2out,minInd,step1=step1,step2=step2)
        if isnan(varxy)
            flag += 2^3 #var off grid
        end
        if (varx<=0) | (vary<=0)
            flag += 2^4 #bad error curvature
        end
        if (varxx<=0) | (varyy<=0)
            flag += 2^5 #very bad error curvature
        end
        return (((minVal_interp[1], minVal_interp[2]), minChi_interp, (valrng1[minInd[1]], valrng2[minInd[2]]), minChi, minInd, flag, varx, vary, varxy, varxx, varyy), (valrng1, valrng2), chi2out)
    end
end

function sampler_2d_hierarchy_var(chi2_fun,lvltup;step1=1,step2=1,minres1=1//10,minres2=1//100)
    lvllen = length(lvltup)
    out = []

    global_minVal_interp1 = 0
    global_minVal_interp2 = 0
    global_minChi_interp = 0
    global_minVal1 = 0
    global_minVal2 = 0
    global_minChi = 0
    global_minInd = CartesianIndex(1,1)
    global_flag = 0
    for (lvlind,lvl) in enumerate(lvltup)
        if lvlind == 1
            uprng1 = lvl[1]
            uprng2 = lvl[2]
        else
            uprng1 = shift_range(lvl[1],round_step(global_minVal_interp1,minres1))
            uprng2 = shift_trim_range(lvl[2],round_step(global_minVal_interp2,minres2))
        end
        (((minVal_interp1, minVal_interp2), minChi_interp, (minVal1, minVal2), minChi, minInd, flag), valrng, chi2out) = sampler_2d_dense(chi2_fun,(uprng1,uprng2))
        push!(out,(((minVal_interp1, minVal_interp2), minChi_interp, (minVal1, minVal2), minChi, minInd, flag), valrng, chi2out))
        if minChi < global_minChi # this should not be the interpolated value, it should be the min measured chi2 value
            (global_minVal_interp1, global_minVal_interp2), global_minChi_interp, (global_minVal1, global_minVal2), global_minChi, global_minInd, global_flag = (minVal_interp1, minVal_interp2), minChi_interp, (minVal1, minVal2), minChi, minInd, flag
        end
        if lvlind == lvllen
            if lvlind == 1
                varx, vary, varxy, varxx, varyy = err2d(lvl,chi2out,minInd;step1=step1,step2=step2)
            else
                varx, vary, varxy, varxx, varyy = err2d((uprng1,uprng2),chi2out,minInd;step1=step1,step2=step2)
            end
            if isnan(varx)
                global_flag += 2^3 #var off grid
            end
            if (varx<=0) | (vary<=0)
                global_flag += 2^4 #bad error curvature
            end
            if (varxx<=0) | (varyy<=0)
                global_flag += 2^5 #very bad error curvature
            end
        return (((global_minVal_interp1, global_minVal_interp2), global_minChi_interp, (global_minVal1, global_minVal2), global_minChi, global_minInd, global_flag, varx, vary, varxy, varxx, varyy), out)
        end
    end
end

function err1d(vrng,chi2,minind;stepx=1)
    lval = length(vrng)
    subgrid = minind .+ (-stepx:stepx)
    if (subgrid[1] < 1) | (lval < subgrid[end])
        return NaN
    else
        c1, c2, c3 = chi2[subgrid]
        x1, x2, x3 = vrng[subgrid]
        dx = (x2-x1)
        if !((x3-x2) ≈ dx)
            @warn "Non uniform grid error estimates not implemented"
        end
        return (dx^2)/(c1-2*c2+c3)
    end
end

# we could convert step to units of the scan... think about
# the right units for that
function err2d(vrngs,chi2,minind;step1=1,step2=1)
    vrng1, vrng2 = vrngs
    lval1 = length(vrng1)
    lval2 = length(vrng2)
    subgrid1 = minind[1] .+ (-step1:step1)
    subgrid2 = minind[2] .+ (-step2:step2)
    # we could consider defaulting to the 1d if one direction is on grid

    if (subgrid1[1] < 1) | (lval1 < subgrid1[end]) | (subgrid2[1] < 1) | (lval2 < subgrid2[end])
        return NaN, NaN, NaN, NaN, NaN
    else
        locgrid = chi2[subgrid1,subgrid2]
        x1, x2, x3 = vrng1[subgrid1]
        dx = (x2-x1)
        if !((x3-x2) ≈ dx)
            @warn "Non uniform grid error estimates not implemented"
        end
        y1, y2, y3 = vrng2[subgrid2]
        dy = (y2-y1)
        if !((y3-y2) ≈ dy)
            @warn "Non uniform grid error estimates not implemented"
        end
        Hxx = (locgrid[3,2]-2*locgrid[2,2]+locgrid[1,2])/(dx^2)
        Hyy = (locgrid[2,3]-2*locgrid[2,2]+locgrid[2,1])/(dy^2)
        Hxy = (locgrid[3,3]-locgrid[3,1]-locgrid[1,3]+locgrid[1,1])/(4*dx*dy)
        det0 = Hxx*Hyy-Hxy^2
        errx= Hyy/det0
        erry= Hxx/det0
        errxy= -Hxy/det0
        return errx, erry, errxy, 1/Hxx, 1/Hyy
    end
end
