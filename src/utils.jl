# Utils Module
using ParallelDataTransfer

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

# Task-Affinity CPU Locking in multinode SlurmContext
function slurm_cpu_lock()
    getinfo_worker(workerid::Int) = @getfrom workerid myid(), ThreadPinning.sched_getcpu(), gethostname()
    idlst = getinfo_worker.(workers()); df = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]); push!(df,idlst...)
    gdf = groupby(df,:hostname)
    @spawnat 1 ThreadPinning.pinthread(1)
    for sgdf in gdf, (sindx, sworker) in enumerate(sgdf.workerid)
        sendto(sworker, sindx=sindx)
        @spawnat sworker ThreadPinning.pinthread(sindx-1)
    end
    # Helpful Worker Info Printing
    idlst = getinfo_worker.(workers()); df = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]); push!(df,idlst...)
    gdf = groupby(df,:hostname); dfc = combine(gdf, nrow, :workerid => minimum, :workerid => maximum, :physcpu => minimum, :physcpu => maximum)
    println("$(gethostname()) running Main on worker: $(myid()) cpu: $(ThreadPinning.sched_getcpu())")
    for row in Tables.namedtupleiterator(dfc)
        println("$(row.hostname) running $(row.nrow) workers: $(row.workerid_minimum)->$(row.workerid_maximum) cpus: $(row.physcpu_minimum)->$(row.physcpu_maximum)")
    end
    flush(stdout)
end
