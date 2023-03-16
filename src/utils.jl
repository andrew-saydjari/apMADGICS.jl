# Utils Module
nanmean(x) = mean(filter(!isnan,x))
nanmean(x,y) = mapslices(nanmean,x,dims=y)

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

nanmedian(x) = median(filter(!isnan,x))
nanmedian(x,y) = mapslices(nanmedian,x,dims=y)

function gaussian_post(x,x0,s0)
    return 1/sqrt(2π)/s0*exp(-0.5*((x-x0)/s0)^2)
end

function sqrt_nan(x)
    if x < 0
        return NaN
    else
        return sqrt(x)
    end
end

function issquare(X)
    sz = size(X)
    return sz[1] == sz[2]
end

function covellipse(Σ; μ=[0,0], n_std=1, n_ellipse_vertices=100)
    λ, U = eigen(Σ)
    S = n_std * U * diagm(.√λ)
    θ = range(0, 2π; length = n_ellipse_vertices)
    A = S * [cos.(θ)'; sin.(θ)']
    return [μ[1] .+ A[1, :], μ[2] .+ A[2, :]]
end

function tuple1dprint(x)
    out = Int[]
    for i=1:length(x)
        push!(out,length(x[i]))
    end
    push!(out,sum(out))
    println(out)
end

function tuple2dprint(x)
    out = Int[]
    for i=1:length(x)
        push!(out,length(x[i][1])*length(x[i][2]))
    end
    push!(out,sum(out))
    println(out)
end

function nanify(x,msk)
    out = zeros(eltype(x),length(msk))
    out[msk] .= x
    out[.!msk] .= NaN
    return out
end
