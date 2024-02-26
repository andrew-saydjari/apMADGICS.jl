# Utils Module
using ParallelDataTransfer, Distributed, Suppressor, SparseArrays

function isnanorzero(x)
    return isnan(x) | iszero(x)
end

nanzeromean(x) = if all(isnanorzero,x)
    NaN
else
    mean(filter(!isnanorzero,x))
end
nanzeromean(x,y) = mapslices(nanzeromean,x,dims=y)

nansum(x) = sum(filter(!isnan,x))
nansum(x,y) = mapslices(nansum,x,dims=y)

nanzerosum(x) = if all(isnanorzero,x)
    NaN
else
    sum(filter(!isnanorzero,x))
end
nanzerosum(x,y) = mapslices(nanzerosum,x,dims=y)

nanzeromedian(x) = if all(isnanorzero,x)
    NaN
else
    median(filter(!isnanorzero,x))
end
nanzeromedian(x,y) = mapslices(nanzeromedian,x,dims=y)

nanzeroiqr(x) = if all(isnanorzero,x)
    NaN
else
    iqr(filter(!isnanorzero,x))/1.34896
end
nanzeroiqr(x,y) = mapslices(nanzeroiqr,x,dims=y)

function inv_nan(A)
    if any(isnan,A)
        return NaN*ones(size(A))
    else
        return inv(A)
    end
end

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

function log10s(x)
    if x<=0
        return NaN
    else
        log10(x)
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
    dfo = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]);
    activateout = @capture_out begin
        idlst = getinfo_worker.(workers()); df = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]); push!(df,idlst...)
        gdf = groupby(df,:hostname)
        @spawnat 1 ThreadPinning.pinthread(0)
        for sgdf in gdf, (sindx, sworker) in enumerate(sgdf.workerid)
            sendto(sworker, sindx=sindx)
            @spawnat sworker ThreadPinning.pinthread(sindx-1)
        end
        rmprocs(2)
        idlst = getinfo_worker.(workers()); push!(dfo,idlst...)
    end
    # Helpful Worker Info Printing
    gdf = groupby(dfo,:hostname); dfc = combine(gdf, nrow, :workerid => minimum, :workerid => maximum, :physcpu => minimum, :physcpu => maximum)
    println("$(gethostname()) running Main on worker: $(myid()) cpu: $(ThreadPinning.sched_getcpu())")
    for row in Tables.namedtupleiterator(dfc)
        println("$(row.hostname) running $(row.nrow) workers: $(row.workerid_minimum)->$(row.workerid_maximum) cpus: $(row.physcpu_minimum)->$(row.physcpu_maximum)")
    end
    flush(stdout)
end

function v2z(v; c=299792.458)
    return sqrt((1+v/c)/(1-v/c))-1
end

function z2v(z; c=299792.458)
    return ((z+1)^2-1)/((z+1)^2+1)*c
end

function z2pix(z; delLog=6e-6)
    return log10(z+1)/delLog
end

function pix2v(x; delLog=6e-6)
    z = 10^(x*delLog)-1
    z2v(z)
end

function prop_p2z(p; delLog=6e-6)
    return log(10)*10^(p*delLog) 
 end
 
 function prop_z2v(z)
    return abs(4*(z+1)/((z+1)^2+1))
 end

 function initalize_git(git_dir)
    git_commit = LibGit2.head(git_dir)
    git_repo = LibGit2.GitRepo(git_dir)
    git_head = LibGit2.head(git_repo)
    git_branch = LibGit2.shortname(git_head)
    println("Running on branch: $git_branch, commit: $git_commit"); flush(stdout)
    return git_branch, git_commit
 end

 function grow_msk2d(msk; rad=1)
    (sx, sy) = size(msk)
    msknew = zeros(Bool,(sx,sy))
    for i=1:sx
        for j=1:sy
            srngx = maximum([1,i-rad]):minimum([sx,i+rad])
            srngy = maximum([1,j-rad]):minimum([sy,j+rad])
            msknew[i,j] = any(msk[srngx,srngy])
        end
    end
    return msknew
end

function chi2red_fluxscale(chi2r, flux; fc=0.0)
    return chi2r/(1 + (fc*flux)^2)
end

function instrument_lsf_sparse_matrix(λ_input, λ_output, R)
    row, col, val = Int[], Int[], Float64[]
    for (indx, λ) in enumerate(λ_output)
        msk, ϕ = instrumental_lsf_kernel(λ_input, λ, R)
        push!(row,(indx.*ones(Int,count(msk)))...)
        push!(col,findall(msk)...)
        push!(val, ϕ...)
    end
    return sparse(row,col,val)
end

function instrumental_lsf_kernel(λ, λc, R)
    σ, (lower, upper) = lsf_sigma_and_bounds(λc, R)
    msk = (lower .<= λ .<= upper)
    ϕ = exp.(-((λ[msk].-λc).^2)/(2σ^2))
    ϕ ./= dropdims(sum(ϕ,dims=1),dims=1)
    return msk, ϕ
end

function lsf_sigma_and_bounds(λ, R; σ_window=10)
     σ = lsf_sigma(λ, R)
    return σ, (λ - σ_window * σ, λ + σ_window * σ)
end

fwhm2sigma = 1/(2 * sqrt(2 * log(2)))
function lsf_sigma(λ, R)
    return (λ / R) * fwhm2sigma
end