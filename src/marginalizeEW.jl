## Marginalization Module (needs to wrap a bit more)

function marginalize_flux_err(chi2lst,fluxlst,dfluxlst,refchi2val;margin_len=361)
    pweight = exp.(0.5*(-chi2lst.+refchi2val))
    pweight ./= sum(filter(.!isnan,pweight))

    exval = vcat(fluxlst.-dfluxlst.*3, fluxlst.+dfluxlst.*3)
    if count(.!isnan.(exval)) > 1
        minx, maxx = extrema(filter(.!isnan,exval))
        #I am not sure how decide on the length of this range obj
        xrng = range(minx,maxx,length=margin_len) 
        
        new_gauss = zeros(margin_len);
        for indlst = 1:length(fluxlst)
            est_flux = fluxlst[indlst]
            sig_flux = dfluxlst[indlst]
            if ((.!isnan(est_flux)) & (.!isnan(sig_flux)))
                new_gauss .+= gaussian_post.(xrng,Ref(est_flux),Ref(sig_flux))*pweight[indlst]
            end
        end
        wvec = ProbabilityWeights(new_gauss)
        out = (mean(xrng,wvec), std(xrng,wvec))
    else
        out = (NaN, NaN)
    end
    return out
end

# are we wasting calls by doing the components and delchi2 seperately?
# can fix this on the next pass...
# should we be doing multiplicative refinement inside this loop?
# need to make deblend and posterior asym component ordering consistent
function sample_chi2_flux_dflux(samp_lst,intup;waveaxis=wavetarg,delLogSpace=delLog)
    lenit = length(samp_lst)
    chi2lst = zeros(lenit)
    fluxlst = zeros(lenit)
    dfluxlst = zeros(lenit)
    (Rs,Ctotinv0,Xd_obs,wave_obs,Dscale,Vcomb_0,V_dib,pre_Vslice,dib_center,scan_offset) = intup
    chi2_wrapper_partial = Base.Fix2(chi2_wrapper2d,(Rs,Ctotinv0,Xd_obs,wave_obs,Dscale,Vcomb_0,V_dib,pre_Vslice,dib_center,scan_offset))
    
    for (samp_ind, samp_tup) in enumerate(samp_lst)
        chi2lst[samp_ind] = chi2_wrapper_partial(samp_tup)
        Ctotinv, Vcomb, V_dibc, V_dibr = update_Ctotinv_Vdib(samp_tup,Ctotinv0.matList[1],Rs,Dscale,Vcomb_0,V_dib,scan_offset)
        
        x_comp_lst = deblend_components_all_asym(Ctotinv, Xd_obs, (V_dibr,), (V_dibc,))
        fluxlst[samp_ind] = sum(x_comp_lst[1].*waveaxis)*delLogSpace
        
        Cii = get_posterior_from_prior_asym(Ctotinv,V_dibc,V_dibr)
        dfluxlst[samp_ind] = sqrt_nan(waveaxis'*(Cii*waveaxis))*delLogSpace
    end
    
    return chi2lst, fluxlst, dfluxlst
end
