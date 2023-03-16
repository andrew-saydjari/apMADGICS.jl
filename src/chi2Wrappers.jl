## Problem specific chi2 wrappers
round_step(x, step) = round(Int,x / step) * step
indInt(sval) = round(Int,fld1(10 *sval +6,10)-1);
indTenth(sval) = mod1(round(Int,mod1(10 *sval +6,10)),10);

sigrng = 4//10:1//100:4//1
sigstep0 = step(sigrng)
minoffset0  = round(Int,sigrng[1] / sigstep0)-1
sigScanFun(x; step = sigstep0, minoffset = minoffset0)  = round(Int,x / step) .-minoffset

# passing Vcomb_0 is new to some of the codes, we need to revisit having to pass these intups
function chi2_wrapper2d(svals,intup;sigslice=4)
    sval1, sval2 = svals
    (simplemsk,Ctotinv,Xd_obs,wave_obs,Dscale,Vcomb_0,V_new,new_center) = intup
    # transform center shift to index
    rval = indInt(sval1)
    tval = indTenth(sval1)
    # transform sigma scan to index
    sigindx = sigScanFun(sval2)
    # calculate window where DIB has support (speed up computation)
    lindx = findfirst(wave_obs .>= (new_center*10^(sval1*6e-6) - sigslice*sval2))
    rindx = findlast(wave_obs .<= (new_center*10^(sval1*6e-6) + sigslice*sval2))
    slrng = lindx:rindx
    return woodbury_update_inv_tst_sub(
        Ctotinv,
        Xd_obs,
        Dscale*(ShiftedArrays.circshift(view(V_new,:,:,sigindx,tval),(rval,0)))[simplemsk,:],
        slrng
        #combine Dscale with Ctotinv when accepting a V_new?
    )
end

function chi2_wrapper(sval,intup)
    (simplemsk,Ctotinv,Xd_obs,Dscale,V_new) = intup
    # transform val to index
    rval = indInt(sval)
    tval = indTenth(sval)
    return woodbury_update_inv_tst(
        Ctotinv,
        Xd_obs,
        Dscale*(ShiftedArrays.circshift(view(V_new,:,:,tval),(rval,0)))[simplemsk,:] 
        #combine Dscale with Ctotinv when accepting a V_new?
    )
end

function update_Ctotinv_Vstarstarlines(svald,Ainv,simplemsk,Dscale,Vcomb_0,V_starlines)
    rval = indInt(svald)
    tval = indTenth(svald)
    V_starlines_c = circshift(view(V_starlines,:,:,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_starlines_r = Dscale*(ShiftedArrays.circshift(view(V_starlines,:,:,tval),(rval,0)))[simplemsk,:]
    Vcomb = hcat(Vcomb_0,V_starlines_r);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_starlines_c, V_starlines_r
end

function update_Ctotinv_Vdib(samp_tup,Ainv,simplemsk,Dscale,Vcomb_0,V_dib)
    (svald,sigvald) = samp_tup
    rval = indInt(svald)
    tval = indTenth(svald)
    sigindx = sigScanFun(sigvald)
    V_dibc = circshift(view(V_dib,:,:,sigindx,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_dibr = Dscale*ShiftedArrays.circshift(view(V_dib,:,:,sigindx,tval),(rval,0))[simplemsk,:]
    Vcomb = hcat(Vcomb_0,V_dibr);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_dibc, V_dibr
end