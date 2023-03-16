## GridSearch Module

function shift_range(inrng,shift)
    start, stop, step1 = inrng[1], inrng[end], step(inrng)
    return range(start+shift,stop+shift,step=step1)
end

function shift_trim_range(inrng,shift; minv=4//10, maxv=4)
    start, stop, step1 = inrng[1], inrng[end], step(inrng)
    return range(maximum((start+shift,minv)),minimum((stop+shift,maxv)),step=step1)
end

function quadratic_interp3(x,y)
    x1x2 = x[1]-x[2]
    x1x3 = x[1]-x[3]
    x2x3 = x[2]-x[3]
    k1 = y[1]/(x1x2*x1x3)
    k2 = y[2]/(-x1x2*x2x3)
    k3 = y[3]/(x2x3*x1x3)
    xopt = k1*(x[2]+x[3])+k2*(x[1]+x[3])+k3*(x[2]+x[1])
    xopt /= (2*(k1+k2+k3))
    yopt = k1*(xopt-x[2])*(xopt-x[3])+k2*(xopt-x[1])*(xopt-x[3])+k3*(xopt-x[2])*(xopt-x[1])
    if yopt > y[2]
        flag = 2^0
        xoptr = x[2]
        yoptr = y[2]
    else
        flag = 0
        xoptr = xopt
        yoptr = yopt
    end
    return xopt, yopt, flag
end

function quadratic_interp3_2d(vrngs,chi2,minind;step1=1,step2=1)
    vrng1, vrng2 = vrngs
    lval1 = length(vrng1)
    lval2 = length(vrng2)
    subgrid1 = minind[1] .+ (-step1:step1)
    subgrid2 = minind[2] .+ (-step2:step2)

    if (lval1 == 1) & (lval2 != 1)
        xopt, yopt, flag = quadratic_interp3(vrng2[subgrid2],chi2[1,subgrid2])
        return (vrng1[minind[1]], xopt), yopt, flag
    elseif (lval1 != 1) & (lval2 == 1)
        xopt, yopt, flag = quadratic_interp3(vrng1[subgrid1],chi2[subgrid1,1])
        return (xopt, vrng2[minind[2]]), yopt, flag
    elseif (lval1 != 1) & (lval2 != 1)
        locgrid = chi2[subgrid1,subgrid2]
        x1, x2, x3 = vrng1[subgrid1]
        dx = (x2-x1)
        if !((x3-x2) ≈ dx)
            warn("Non uniform grid error estimates not implemented")
        end
        y1, y2, y3 = vrng2[subgrid2]
        dy = (y2-y1)
        if !((y3-y2) ≈ dy)
            warn("Non uniform grid error estimates not implemented")
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

        minsval =  x0 .- (Hinv*G)
        function f(x)
            dx = (x-x0)
            return locgrid[2,2] + G'*dx + (dx'*(H*dx))/2
        end
        minvali = f(minsval)

        # In 1d, this is definitely not possible
        # in 2d, this flag might be important
        if minvali > locgrid[2,2]
            flag = 2^0
            minopt = minind
            valopt = locgrid[2,2]
        else
            flag = 0
            minopt = minsval
            valopt = minvali
        end

        return minopt, valopt, flag
    else
        warn("I can't interpolate a 0d scan...")
    end
end

function sampler_1d_dense(chi2_fun,valrng)
    lval = length(valrng)
    chi2out = zeros(lval)
    for (ind,val) in enumerate(valrng)
        chi2out[ind] = chi2_fun(val)
    end
    minval, minind = findmin(chi2out)
    subgrid = minind .+ (-1:1)
    if (minind == 1) | (minind == lval)
        flag = 2^1 # edge of grid bit
        return ((minval, valrng[minind], valrng[minind], minind, flag), valrng, chi2out)
    else
        minvali, minsvali, flag = quadratic_interp3(valrng[subgrid],chi2out[subgrid])
        #need to improve the notation with "s" and no "s" on the mins
        return ((minsvali, minvali, valrng[minind], minind, flag), valrng, chi2out)
    end
end

function sampler_1d_dense_var(chi2_fun,valrng;stepx=1)
    chi2out = zeros(length(valrng))
    for (ind,val) in enumerate(valrng)
        chi2out[ind] = chi2_fun(val)
    end
    minval, minind = findmin(chi2out)
    subgrid = minind .+ (-1:1)
    minvali, minsvali, flag = quadratic_interp3(valrng[subgrid],chi2out[subgrid])
    varx = err1d(valrng,chi2out,minind;stepx=stepx)
    if isnan(varx)
        flag += 2^2 #var off grid
    end
    if varx<=0
        flag += 2^3 #bad error curvature
    end
    return ((minsvali, minvali, valrng[minind], minind, flag, varx), valrng, chi2out)
end

function sampler_1d_hierarchy_var(chi2_fun,lvltup;minres=1//10,stepx=1)
    lvllen = length(lvltup)
    out = []

    globalminvali = 0
    globalminsvali = NaN
    gloablminsval = NaN
    globalminind = 1
    globalflag = 0
    for (lvlind,lvl) in enumerate(lvltup)
        if lvlind == 1
            ((minvali, minsvali, minsval, minind, flag), valrng, chi2out) = sampler_1d_dense(chi2_fun,lvl)
            push!(out,((minvali, minsvali, minsval, minind, flag), valrng, chi2out))
            globalminvali, globalminsvali, gloablminsval, globalminind, globalflag = minvali, minsvali, minsval, minind, flag
        else
            uprng = shift_range(lvl,round_step(globalminsvali,minres))
            ((minvali, minsvali, minsval, minind, flag), valrng, chi2out)  = sampler_1d_dense(chi2_fun,uprng)
            push!(out,((minvali, minsvali, minsval, minind, flag), valrng, chi2out))
            if minvali < globalminvali # this was wrong and svali before... that is concerning...
                globalminvali, globalminsvali, gloablminsval, globalminind, globalflag = minvali, minsvali, minsval, minind, flag
            end
        end
        if lvlind == lvllen
            varx = err1d(lvl,chi2out,minind,stepx=stepx)
                if isnan(varx)
                    globalflag += 2^2 #var off grid
                end
                if varx<=0
                    globalflag += 2^3 #bad error curvature
                end
            return ((globalminvali, globalminsvali, gloablminsval, globalminind, globalflag, varx), out)
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
    minval, minind = findmin(chi2out)
    # now we have a gridx and gridy edge bit (might need to update bit identities)
    flag = (2^1)*((minind[1] == 1) | (minind[1] == lval1))*(lval1!=1) + (2^2)*((minind[2] == 1) | (minind[2] == lval2))*(lval2!=1)
    if flag !=0
        return ((minval, (valrng1[minind[1]], valrng2[minind[2]]), (valrng1[minind[1]], valrng2[minind[2]]), minind, flag), valrngtup, chi2out)
    else
        minsvali, minvali, flag = quadratic_interp3_2d(valrngtup,chi2out,minind;step1=1,step2=1)
        return ((minvali, (minsvali[1], minsvali[2]), (valrng1[minind[1]], valrng2[minind[2]]), minind, flag), valrngtup, chi2out)
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
    minval, minind = findmin(chi2out)
    # now we have a gridx and gridy edge bit (might need to update bit identities)
    flag = (2^1)*((minind[1] == 1) | (minind[1] == lval1)) + (2^2)*((minind[2] == 1) | (minind[2] == lval2))
    if flag !=0
        return ((minval, (valrng1[minind[1]], valrng2[minind[2]]), (valrng1[minind[1]], valrng2[minind[2]]), minind, flag), valrngtup, chi2out)
    else
        minsvali, minvali, flag = quadratic_interp3_2d(vrngs,chi2,minind;step1=1,step2=1)
        return ((minvali, (minsvali[1], minsvali[2]), (valrng1[minind[1]], valrng2[minind[2]]), minind, flag), valrngtup, chi2out)
    end
    # we may want 1d quadratic interpolations here on if cases

    # one could do 1d error bars on edge of grid problems with careful if cases
    varx, vary, varxy, varxx, varyy = err2d((valrng1,valrng2),chi2out,minind,step1=step1,step2=step2)
    if isnan(varxy)
        flag += 2^3 #var off grid
    end
    if (varx<=0) | (vary<=0)
        flag += 2^4 #bad error curvature
    end
    if (varxx<=0) | (varyy<=0)
        flag += 2^5 #very bad error curvature
    end
    return ((minval, (valrng1[minind[1]], valrng2[minind[2]]), (valrng1[minind[1]], valrng2[minind[2]]), minind, flag, varx, vary, varxy, varxx, varyy), valrng1, valrng2, chi2out)
end

function sampler_2d_hierarchy_var(chi2_fun,lvltup;step1=1,step2=1,minres1=1//10,minres2=1//100)
    lvllen = length(lvltup)
    out = []

    globalminvali = 0
    globalminsvali1 = NaN
    globalminsvali2 = NaN
    gloablminsval1 = NaN
    gloablminsval2 = NaN
    globalminind = CartesianIndex(1,1)
    globalflag = 0
    for (lvlind,lvl) in enumerate(lvltup)
        if lvlind == 1
            ((minvali, (minsvali1, minsvali2), (minsval1, minsval2), minind, flag), valrng, chi2out) = sampler_2d_dense(chi2_fun,lvl)
            push!(out,((minvali, (minsvali1, minsvali2), (minsval1, minsval2), minind, flag), valrng, chi2out))
            globalminvali, (globalminsvali1, globalminsvali2), (gloablminsval1, gloablminsval2), globalminind, globalflag = minvali, (minsvali1, minsvali2), (minsval1, minsval2), minind, flag
        else
            uprng1 = shift_range(lvl[1],round_step(globalminsvali1,minres1))
            uprng2 = shift_trim_range(lvl[2],round_step(globalminsvali2,minres2))
            ((minvali, (minsvali1, minsvali2), (minsval1, minsval2), minind, flag), valrng, chi2out) = sampler_2d_dense(chi2_fun,(uprng1,uprng2))
            push!(out,((minvali, (minsvali1, minsvali2), (minsval1, minsval2), minind, flag), valrng, chi2out))
            if minvali < globalminvali
                globalminvali, (globalminsvali1, globalminsvali2), (gloablminsval1, gloablminsval2), globalminind, globalflag = minvali, (minsvali1, minsvali2), (minsval1, minsval2), minind, flag
            end
        end
        if lvlind == lvllen
            if lvlind == 1
                varx, vary, varxy, varxx, varyy = err2d(lvl,chi2out,minind;step1=step1,step2=step2)
            else
                varx, vary, varxy, varxx, varyy = err2d((uprng1,uprng2),chi2out,minind;step1=step1,step2=step2)
            end
            if isnan(varx)
                globalflag += 2^3 #var off grid
            end
            if (varx<=0) | (vary<=0)
                globalflag += 2^4 #bad error curvature
            end
            if (varxx<=0) | (varyy<=0)
                globalflag += 2^5 #very bad error curvature
            end
        return ((globalminvali, (globalminsvali1, globalminsvali2), (gloablminsval1, gloablminsval2), globalminind, globalflag, varx, vary, varxy, varxx, varyy), out)
        end
    end
end

function err1d(vrng,chi2,minind;stepx=1)
    lval = length(vrng)
    subgrid = minind .+ (-stepx:stepx)
    if (subgrid[1] < 1) | (lval < subgrid[end])
        return NaN
    else
        chi2_1, chi2_2, chi2_3 = chi2[subgrid]
        x1, x2, x3 = vrng[subgrid]
        dx = (x2-x1)
        if !((x3-x2) ≈ dx)
            warn("Non uniform grid error estimates not implemented")
        end
        return (dx^2)/(chi2_1-2*chi2_2+chi2_3)
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
            warn("Non uniform grid error estimates not implemented")
        end
        y1, y2, y3 = vrng2[subgrid2]
        dy = (y2-y1)
        if !((y3-y2) ≈ dy)
            warn("Non uniform grid error estimates not implemented")
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
