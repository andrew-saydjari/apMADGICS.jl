## Problem specific chi2 wrappers
round_step(x, step) = round(Int,x / step) * step
indInt(sval) = round(Int,fld1(10 *sval + 55//10,10)-1);
indTenth(sval) = mod1(round(Int,mod1(10 *sval +6,10)),10);

sigrng = 4//10:1//100:4//1
sigstep0 = step(sigrng)
minoffset0  = round(Int,sigrng[1] / sigstep0)-1
sigScanFun(x; step = sigstep0, minoffset = minoffset0)  = round(Int,x / step) .-minoffset

#MDispatch to handle nothing finds (suspect missing chip case)
function slicer(lindx::Int,rindx::Int,widx)
    return lindx:rindx
end
function slicer(lindx::Int,rindx::Nothing,widx)
    return lindx:(lindx+widx)
end
function slicer(lindx::Nothing,rindx::Int,widx)
    return (rindx-widx):rindx
end

# passing Vcomb_0 is new to some of the codes, we need to revisit having to pass these intups
function chi2_wrapper2d(svals,intup;sigslice=4)
    sval1, sval2 = svals
    (simplemsk,Ctotinv,Xd_obs,wave_obs,Dscale,Vcomb_0,V_new,pre_Vslice,new_center,center_offset) = intup
    # transform center shift to index
    rval = indInt(sval1+center_offset)
    tval = indTenth(sval1+center_offset)
    # transform sigma scan to index
    sigindx = sigScanFun(sval2)
    # calculate window where DIB has support (speed up computation)
    lindx = findfirst(wave_obs .>= (new_center*10^(sval1*6e-6) - sigslice*sval2))
    rindx = findlast(wave_obs .<= (new_center*10^(sval1*6e-6) + sigslice*sval2))
    widx = round(Int,sigslice*sval2/0.22) # this is a hardcoded approximate pixel size in Angstroms
    slrng = slicer(lindx,rindx,widx)
    pre_Vslice .= view(ShiftedArrays.circshift(view(V_new,:,:,sigindx,tval),(rval,0)),simplemsk,:)
    pre_Vslice .*= Dscale 
    return woodbury_update_inv_tst_sub(
        Ctotinv,
        Xd_obs,
        pre_Vslice,
        slrng
    )
end

function chi2_wrapper(sval,intup)
    (simplemsk,Ctotinv,Xd_obs,Dscale,V_new,pre_Vslice) = intup
    # transform val to index
    rval = indInt(sval)
    tval = indTenth(sval)
    pre_Vslice .= view(ShiftedArrays.circshift(view(V_new,:,:,tval),(rval,0)),simplemsk,:)
    pre_Vslice .*= Dscale 
    return woodbury_update_inv_tst(
        Ctotinv,
        Xd_obs,
        pre_Vslice
    )
end

function update_Ctotinv_Vstarstarlines(svald,Ainv,simplemsk,Dscale,Vcomb_0,V_starlines)
    rval = indInt(svald)
    tval = indTenth(svald)
    V_starlines_c = circshift(view(V_starlines,:,:,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_starlines_r = ShiftedArrays.circshift(view(V_starlines,:,:,tval),(rval,0))[simplemsk,:]
    V_starlines_r .*= Dscale
    Vcomb = hcat(Vcomb_0,V_starlines_r);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_starlines_c, V_starlines_r
end

function update_Ctotinv_Vstarstarlines_asym(svald,Ainv,simplemsk,Dscale,Vcomb_0,V_starlines,V_starlines_refLSF)
    rval = indInt(svald)
    tval = indTenth(svald)
    V_starlines_c = circshift(view(V_starlines_refLSF,:,:,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_starlines_r = ShiftedArrays.circshift(view(V_starlines,:,:,tval),(rval,0))[simplemsk,:]
    V_starlines_r .*= Dscale
    Vcomb = hcat(Vcomb_0,V_starlines_r);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_starlines_c, V_starlines_r
end

function update_Ctotinv_Vstarstarlines2_asym(svald,Ainv,simplemsk,Dscale,Vcomb_0,V_starlines,V_starlines_refLSF,V_cor)
    rval = indInt(svald)
    tval = indTenth(svald)

    V_starlines_c = circshift(view(V_starlines_refLSF,:,:,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_starlines_r = ShiftedArrays.circshift(view(V_starlines,:,:,tval),(rval,0))[simplemsk,:]
    V_starlines_r .*= Dscale

    V_cor_c = circshift(view(V_cor,:,:,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_cor_r = ShiftedArrays.circshift(view(V_cor,:,:,tval),(rval,0))[simplemsk,:]
    V_cor_r .*= Dscale

    Vcomb = hcat(Vcomb_0,V_starlines_r,V_cor_r);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_starlines_c, V_starlines_r, V_cor_c, V_cor_r
end

function update_Ctotinv_Vdib(samp_tup,Ainv,simplemsk,Dscale,Vcomb_0,V_dib,scan_offset)
    (svald,sigvald) = samp_tup
    rval = indInt(svald+scan_offset)
    tval = indTenth(svald+scan_offset)
    sigindx = sigScanFun(sigvald)
    V_dibc = circshift(view(V_dib,:,:,sigindx,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_dibr = ShiftedArrays.circshift(view(V_dib,:,:,sigindx,tval),(rval,0))[simplemsk,:]
    V_dibr .*= Dscale
    Vcomb = hcat(Vcomb_0,V_dibr);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_dibc, V_dibr
end

function update_Ctotinv_Vdib_asym(samp_tup,Ainv,simplemsk,Dscale,Vcomb_0,V_dib,V_dib_noLSF,scan_offset)
    (svald,sigvald) = samp_tup
    rval = indInt(svald+scan_offset)
    tval = indTenth(svald+scan_offset)
    sigindx = sigScanFun(sigvald)
    V_dibc = circshift(view(V_dib_noLSF,:,:,sigindx,tval),(rval,0)) # this needs NaNs or something... VBAD
    V_dibr = ShiftedArrays.circshift(view(V_dib,:,:,sigindx,tval),(rval,0))[simplemsk,:]
    V_dibr .*= Dscale
    Vcomb = hcat(Vcomb_0,V_dibr);
    Ctotinv = LowRankMultMat([Ainv,Vcomb],wood_precomp_mult,wood_fxn_mult);
    return Ctotinv, Vcomb, V_dibc, V_dibr
end

