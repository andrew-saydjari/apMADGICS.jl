## Component and Posterior Module

function deblend_components_all(Σ_inv, Xd, prior_lst)
    out = []
    lprior = length(prior_lst)
    
    Σinvx = (Σ_inv*Xd)
    
    for i=1:lprior
        push!(out,get_component_from_prior(prior_lst[i],Σinvx))
    end
    return out
end

function deblend_components_all_asym(Σ_inv, Xd, prior_lst_R, prior_lst_L)
    out = []
    lprior = length(prior_lst_R)
    @assert lprior == length(prior_lst_L)
    Σinvx = (Σ_inv*Xd)
    
    for i=1:lprior
        push!(out,get_component_from_prior_asym(prior_lst_R[i],prior_lst_L[i],Σinvx))
    end
    return out
end
        
function deblend_components_all_tot(Σ_inv, Xd, prior_lst)
    out = []
    lprior = length(prior_lst)
    
    Σinvx = (Σ_inv*Xd)
    
    for i=1:lprior
        push!(out,get_component_from_prior(prior_lst[i],Σinvx))
    end
    push!(out,Xd'*Σinvx)
    return out
end

function deblend_components_all_asym_tot(Σ_inv, Xd, prior_lst_R, prior_lst_L)
    out = []
    lprior = length(prior_lst_R)
    @assert lprior == length(prior_lst_L)
    Σinvx = (Σ_inv*Xd)
    
    for i=1:lprior
        push!(out,get_component_from_prior_asym(prior_lst_R[i],prior_lst_L[i],Σinvx))
    end
    push!(out,Xd'*Σinvx)
    return out
end

function get_component_from_prior(pobj,Σinvx)
    if typeof(pobj) <: Diagonal
        return pobj*Σinvx
    elseif typeof(pobj) <: Vector
        return Diagonal(pobj)*Σinvx
    elseif !issquare(pobj)
        return pobj*(pobj'*Σinvx)
    else
        return pobj*Σinvx
    end
end

# so only active when using a low rank decomp, you can do whatever
# for the L argument in the other cases, though standard should be to replicated
function get_component_from_prior_asym(pobjR,pobjL,Σinvx)
    if typeof(pobjR) <: Diagonal
        return pobjR*Σinvx
    elseif typeof(pobjR) <: Vector
        return Diagonal(pobjR)*Σinvx
    elseif !issquare(pobjR)
        return pobjL*(pobjR'*Σinvx)
    else
        return pobjR*Σinvx
    end
end

function compute_posterior_ij(Σ_inv, prior_lst_i, prior_lst_j)
    out = []
    lprior = length(prior_lst_i)
    
    for i=1:lprior
        push!(out,get_posterior_ij_from_prior(Σ_inv, prior_lst_i[i],prior_lst_j[i]))
    end
    return out
end

function get_posterior_ij_from_prior(Σ_inv, pobj_i, pobj_j)
    if (typeof(pobj_i) <: Diagonal) & (typeof(pobj_j) <: Diagonal)
        return LowRankMultMat([Σ_inv,pobj_i, pobj_j],[],Cij_fxn_dmat_mult);
    elseif (!issquare(pobj_i)) & (typeof(pobj_j) <: Diagonal)
        return LowRankMultMat([Σ_inv,pobj_i, pobj_j],[],Cij_fxn_comp_dmat_mult);
    elseif (typeof(pobj_i) <: Diagonal) & (!issquare(pobj_j))
        return LowRankMultMat([Σ_inv,pobj_i, pobj_j],[],Cij_fxn_dmat_comp_mult);
    elseif (!issquare(pobj_i)) & (!issquare(pobj_j))
        return LowRankMultMat([Σ_inv,pobj_i, pobj_j],[],Cij_fxn_mult);
    else # how many times would you have to call this for this to be a bad idea for a dense matrix n/3?
        return LowRankMultMat([Σ_inv,pobj_i, pobj_j],[],Cij_fxn_dmat_mult);
    end
end

function compute_posteriors(Σ_inv, prior_lst)
    out = []
    lprior = length(prior_lst)
    
    for i=1:lprior
        push!(out,get_posterior_from_prior(Σ_inv, prior_lst[i]))
    end
    return out
end

# can generalize with multiple dispatch
function compute_posteriors_asym(Σ_inv, prior_lst_L, prior_lst_R)
    out = []
    lprior = length(prior_lst_L)
    
    for i=1:lprior
        push!(out,get_posterior_from_prior_asym(Σ_inv, prior_lst_L[i], prior_lst_R[i]))
    end
    return out
end

function get_posterior_from_prior(Σ_inv, pobj)
    if typeof(pobj) <: Diagonal
        return LowRankMultMat([Σ_inv,pobj],Cii_precomp_mult,Cii_fxn_dmat_mult);
    elseif !issquare(pobj)
        return LowRankMultMat([Σ_inv,pobj],[],Cii_fxn_mult);
    else # how many times would you have to call this for this to be a bad idea for a dense matrix n/3?
        return LowRankMultMat([Σ_inv,pobj],[],Cii_fxn_dmat_mult);
    end
end

# only having it do something for the V case... should be fine for now
function get_posterior_from_prior_asym(Σ_inv, pobjL, pobjR)
    if typeof(pobjL) <: Diagonal
        return LowRankMultMat([Σ_inv,pobjL],Cii_precomp_mult,Cii_fxn_dmat_mult);
    elseif !issquare(pobjL)
        return LowRankMultMat([Σ_inv,pobjL, pobjR],[],Cii_fxn_mult_asym);
    else # how many times would you have to call this for this to be a bad idea for a dense matrix n/3?
        return LowRankMultMat([Σ_inv,pobjL],[],Cii_fxn_dmat_mult);
    end
end
                            
function get_diag_posterior_from_prior_asym(Σ_inv, pobjL, pobjR)
    tempDiagMat = LowRankDiagMat([Σ_inv,pobjL,pobjR],Cii_precomp_asym_diag,Cii_diag_asym_map);
    return diag(tempDiagMat)           
end
