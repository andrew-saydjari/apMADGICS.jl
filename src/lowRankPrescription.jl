## LowRankMultMatIP
using LinearAlgebra
import Base: *
import LinearAlgebra: mul!

struct LowRankMultMatIP{T,pT<:Vector{Matrix{Float64}},F<:Function,G<:Function}
    matList::T
    precompList::pT
    multFunc::F
    multFuncIP::G
end

mul!(W::LowRankMultMatIP,B::Union{Vector,AbstractMatrix}) = W.multFuncIP(W.matList,W.precompList,B)
*(W::LowRankMultMatIP,B::Union{Vector,AbstractMatrix}) = W.multFunc(W.matList,W.precompList,B)
(W::LowRankMultMatIP)(B::Union{Vector,AbstractMatrix}) = W.multFunc(W.matList,W.precompList,B)

function prealloc_mat_gen(sAinv,sVold,sVnew)
    AinvVIinvtX = zeros(sVold[2],sVnew[2])
    VAinvVIinvtX = zeros(sVnew)
    out = zeros(sVnew)
    return AinvVIinvtX, VAinvVIinvtX, out
end

function wood_fxn_mult_mat!(matList::Vector{AbstractMatrix{Float64}},
        precompList::Vector{Matrix{Float64}},
        x::Union{Vector{Float64},AbstractMatrix{Float64}}
    )
    Ainv = matList[1]
    V = matList[2]
    AinvVIinv = precompList[1]
    AinvVIinvtX = precompList[2]
    VAinvVIinvtX = precompList[3]
    out = precompList[4]
    
    mul!(AinvVIinvtX,AinvVIinv',x)
    mul!(VAinvVIinvtX,V,AinvVIinvtX,-1,0)
    VAinvVIinvtX .+= x
    mul!(out,Ainv,VAinvVIinvtX)
    return 
end

function wood_fxn_mult(matList::Vector{AbstractMatrix{Float64}},
        precompList::Vector{Matrix{Float64}},
        x::Union{Vector{Float64},AbstractMatrix{Float64}})
    Ainv = matList[1]
    V = matList[2]
    arg1 = precompList[1]
    return Ainv*(x - V*(arg1'*x))
end

function wood_precomp_mult_mat(matList::Vector{AbstractMatrix{Float64}}, sVnew::Tuple{Int,Int})
    Ainv = matList[1]
    V = matList[2]
    AinvV = Ainv*V
    prealloc = prealloc_mat_gen(size(Ainv),size(V),sVnew)
    return [(AinvV)*inv(I+V'*(AinvV)),prealloc[1],prealloc[2],prealloc[3]]
end

## Low Rank Mult Partners Module  

function densify(x)
    len = size(x.matList[end],1)
    tst = zeros(len)
    out = zeros(len,len)
    for i=1:len
        fill!(tst,0)
        tst[i] = 1
        out[i,:] .= x*tst
    end
    return out
end

function woodbury_update_inv_tst(Ainv::LowRankMultMatIP,Xd,V)
    mul!(Ainv,V)
    AinvV = Ainv.precompList[end]
    XdAinvV = reshape(Xd,1,:)*AinvV
    M = V'*AinvV
    dind = diagind(M)
    M[dind] .+= 1
    Minvc = cholesky!(Symmetric(M))
    return -(XdAinvV*(Minvc\XdAinvV'))[1]/2
end

function woodbury_update_inv_split_tst(Ainv::LowRankMultMatIP,Xd,V2,AinvV1,XdAinvV1,V1TAinvV1;chi2offset=0.0)
    mul!(Ainv,V2)
    AinvV2 = Ainv.precompList[end]
    XdAinvV2 = reshape(Xd,1,:)*AinvV2
    V2TAinvV2 = V2'*AinvV2
    V2TAinvV1 = V2'*AinvV1
    M = [V1TAinvV1 V2TAinvV1';V2TAinvV1  V2TAinvV2] # how efficient is this block allocation?
    dind = diagind(M)
    M[dind] .+= 1
    Minvc = cholesky!(Symmetric(M))
    XdAinvV1V2 = hcat(XdAinvV1,XdAinvV2)
    return -(XdAinvV1V2*(Minvc\XdAinvV1V2'))[1]/2-chi2offset
end

function woodbury_update_inv_tst_res(Ainv::LowRankMultMatIP,Xd,V,Cres::Diagonal)
    mul!(Ainv,V)
    AinvV = Ainv.precompList[end]
    XdAinvV = reshape(Xd,1,:)*AinvV
    XdAinvCresAinvV = (Ainv*reshape(Xd,:,1))'*(Cres*AinvV) # super not ideal
    M = V'*AinvV
    dind = diagind(M)
    M[dind] .+= 1
    Minvc = cholesky!(Symmetric(M))
    return -(XdAinvCresAinvV*(Minvc\XdAinvV'))[1]/2
end

function woodbury_update_inv_tst(Ainv::Union{LowRankMultMat,AbstractMatrix},Xd,V)
    AinvV = Ainv*V
    XdAinvV = reshape(Xd,1,:)*AinvV
    M = V'*AinvV
    dind = diagind(M)
    M[dind] .+= 1
    Minvc = cholesky!(Symmetric(M))
    return -(XdAinvV*(Minvc\XdAinvV'))[1]/2
end

function woodbury_update_inv_tst_sub(Ainv,Xd,V,slrng)
    Ainvsub = LowRankMultMat(Ainv.matList,(Ainv.precompList[1],slrng),wood_fxn_mult_sub)
    AinvV = Ainvsub*view(V,slrng,:)
    XdAinvV = reshape(Xd,1,:)*AinvV
    M = V'*AinvV
    dind = diagind(M)
    M[dind] .+= 1
    Minvc = cholesky!(Symmetric(M))
    return -(XdAinvV*(Minvc\XdAinvV'))[1]/2
end

function wood_precomp_mult(matList)
    Ainv = matList[1]
    V = matList[2]
    AinvV = Ainv*V
    if any(.!isfinite.(AinvV))
        return [AinvV]
    else
        return [(AinvV)*inv(I+V'*(AinvV))]
    end
end

function wood_fxn_mult(matList,precompList,x)
    Ainv = matList[1]
    V = matList[2]
    arg1 = precompList[1]
    return Ainv*(x - V*(arg1'*x))
end

function wood_fxn_mult_sub(matList,precompList,x)
    Ainv = matList[1]
    V = matList[2]
    arg1 = precompList[1]
    slrng = precompList[2]
    out = -(Ainv*(V*(view(arg1',:,slrng)*x)))
    out[slrng,:] += diag(Ainv)[slrng].*x
    return out
end

function Cii_fxn_mult(matList,precompList,x)
    Ctotinv = matList[1]
    Vi = matList[2]
    VVtx = Vi*(Vi'*x)
    return VVtx - Vi*(Vi'*(Ctotinv*VVtx))
end

function Cii_fxn_dmat_mult(matList,precompList,x)
    Ctotinv = matList[1]
    C = matList[2]
    Cx = C*x
    return Cx - C*(Ctotinv*Cx)
end

function Cii_fxn_mult_asym(matList,precompList,x)
    Ctotinv = matList[1]
    Vi = matList[2]
    Vir = matList[3]
    Vtx = (Vi'*x)
    return Vi*Vtx - Vi*(Vir'*(Ctotinv*(Vir*Vtx)))
end

function Cij_fxn_mult(matList,precompList,x)
    Ctotinv = matList[1]
    Vi = matList[2]
    Vj = matList[3]
    return -Vi*(Vi'*(Ctotinv*(Vj*(Vj'*x))))
end

function Cij_fxn_dmat_mult(matList,precompList,x)
    Ctotinv = matList[1]
    Ci = matList[2]
    Cj = matList[3]
    return -Ci'*(Ctotinv*(Cj*x))
end

function Cij_fxn_dmat_comp_mult(matList,precompList,x)
    Ctotinv = matList[1]
    Ci = matList[2]
    Vj = matList[3]
    return -Ci'*(Ctotinv*(Vj*(Vj'*x)))
end

function Cij_fxn_comp_dmat_mult(matList,precompList,x)
    Ctotinv = matList[1]
    Vi = matList[2]
    Cj = matList[3]
    return -Vi*(Vi'*(Ctotinv*(Cj*x)))
end
                                                                        
function Cii_precomp_diag(matList)
    Ctotinv = matList[1]
    Vi = matList[2]
    return [Vi'*(Ctotinv*Vi)]
end

function Cii_diag_map(matList,precompList)
    Vi = matList[2]
    arg1 = precompList[1]
    return dropdims(sum(Vi.^2,dims=2),dims=2).-dropdims(sum(Vi'.*(arg1*Vi'),dims=1),dims=1)
end
                                                                                
function Cii_precomp_asym_diag(matList)
    Ctotinv = matList[1]
    Vir = matList[3]
    return [Vir'*(Ctotinv*Vir)]
end

function Cii_diag_asym_map(matList,precompList)
    Vi = matList[2]
    arg1 = precompList[1]
    return dropdims(sum(Vi.^2,dims=2),dims=2).-dropdims(sum(Vi'.*(arg1*Vi'),dims=1),dims=1)
end                                                                   
