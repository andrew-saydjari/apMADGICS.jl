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

function wood_precomp_mult(matList)
    Ainv = matList[1]
    V = matList[2]
    AinvV = Ainv*V
    return [(AinvV)*inv(I+V'*(AinvV))]
end
