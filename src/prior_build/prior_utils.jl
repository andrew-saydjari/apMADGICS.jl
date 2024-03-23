### utility functions for prior building 

# requires medframes in global scope, that is kinda gross AKS (by AKS)
function generate_poly_prior(adjfibindx,medframes)
    fiber = if adjfibindx>300
        adjfibindx-300
    else
        adjfibindx
    end
    x_model = 15000:0.01:17000
    chipc = 2048:-1:1
    chipb = (2048:-1:1).+2048
    chipa = (2048:-1:1).+2*2048
    subchipa = medframes[chipa,fiber];
    mska = (minimum(subchipa) .<= wavetarg .<= maximum(subchipa))
    mskae = (maximum(subchipa) .<= wavetarg)
    subchipb = medframes[chipb,fiber];
    mskb = (minimum(subchipb) .<= wavetarg .<= maximum(subchipb))
    subchipc = medframes[chipc,fiber];
    mskc = (minimum(subchipc) .<= wavetarg .<= maximum(subchipc));
    mskce = (wavetarg .<= minimum(subchipc));

    gapab = (maximum(subchipb) .< wavetarg .< minimum(subchipa))
    gapab_ind = findall(gapab);

    gapbc = (maximum(subchipc) .< wavetarg .< minimum(subchipb))
    gapbc_ind = findall(gapbc);

    linvec_a = zeros(length(mska))
    linvec_a .= mska #.| mskae;
    # linvec_a[gapab_ind] .= range(0,1,length(gapab_ind));

    linvec_b = zeros(length(mskb))
    linvec_b .= mskb
    # linvec_b[gapab_ind] .= range(1,0,length(gapab_ind))
    # linvec_b[gapbc_ind] .= range(0,1,length(gapbc_ind));

    linvec_c = zeros(length(mska))
    linvec_c .= mskc; #.| mskce;
    # linvec_c[gapbc_ind] .= range(1,0,length(gapbc_ind));

    msknall = .!(mska .| mskb .| mskc);

    nbasis = 1
    msklst = [findall(mska),findall(mskb),findall(mskc)]
    cbasis = BasisFunctions.CosineSeries(nbasis)
    Vout = zeros(length(wavetarg),3,nbasis)
    for chipind = 1:3
        mski = msklst[chipind]
        cntmski = length(mski)
        svec0 = range(0,1,length=cntmski)
        for (ind, sval) in enumerate(svec0)
            Vout[mski[ind],chipind,:] .= cbasis(sval)
        end
    end
    nbasis=2
    cbasis = BasisFunctions.ChebyshevT(nbasis)
    Vout1 = zeros(length(wavetarg),nbasis)
    svec1 = range(-1,1,length=length(wavetarg))
    for (ind, sval) in enumerate(svec1)
        Vout1[ind,:] .= cbasis(sval)
    end
    Vout1[msknall,:].=0;

    nbasis=5
    cbasis = BasisFunctions.CosineSeries(nbasis)
    Vout2 = zeros(length(wavetarg),nbasis)
    svec2 = range(0,1,length=length(wavetarg))
    for (ind, sval) in enumerate(svec2)
        Vout2[ind,:] .= cbasis(sval)
    end
    Vout2[msknall,:].=0;
    
    # Linear Pieces
    mskb_start = (minimum(subchipb) .<= wavetarg .<= 15875)
    mskb_start_i = findall(mskb_start)
    nbasis=2
    cbasis = BasisFunctions.ChebyshevT(nbasis)
    Vout3 = zeros(length(wavetarg),nbasis)
    cntmski3 = length(mskb_start_i)
    svec3 = range(1,0,length=cntmski3)
    for (ind, sval) in enumerate(svec3)
        Vout3[mskb_start_i[ind],:] .= cbasis(sval)
    end
    Vout3[msknall,:].=0;
    
    # Linear Pieces
    mskc_end = (15700 .<= wavetarg .<= maximum(subchipc))
    mskc_end_i = findall(mskc_end)
    nbasis=2
    cbasis = BasisFunctions.ChebyshevT(nbasis)
    Vout4 = zeros(length(wavetarg),nbasis)
    cntmski4 = length(mskc_end_i)
    svec4 = range(0,1,length=cntmski4)
    for (ind, sval) in enumerate(svec4)
        Vout4[mskc_end_i[ind],:] .= cbasis(sval)
    end
    Vout4[msknall,:].=0;
    
    # Linear Pieces
    mska_start = (minimum(subchipa) .<= wavetarg .<= 16550)
    mska_start_i = findall(mska_start)
    nbasis=2
    cbasis = BasisFunctions.ChebyshevT(nbasis)
    Vout5 = zeros(length(wavetarg),nbasis)
    cntmski5 = length(mska_start_i)
    svec5 = range(1,0,length=cntmski5)
    for (ind, sval) in enumerate(svec5)
        Vout5[mska_start_i[ind],:] .= cbasis(sval)
    end
    Vout5[msknall,:].=0;
    
    # Linear Pieces
    mska_end = (16910 .<= wavetarg .<= maximum(subchipa))
    mska_end_i = findall(mska_end)
    nbasis=2
    cbasis = BasisFunctions.ChebyshevT(nbasis)
    Vout7 = zeros(length(wavetarg),nbasis)
    cntmski7 = length(mska_end_i)
    svec7 = range(0,1,length=cntmski7)
    for (ind, sval) in enumerate(svec7)
        Vout7[mska_end_i[ind],:] .= cbasis(sval)
    end
    Vout7[msknall,:].=0;
    
    nbasis = 4
    cbasis = BasisFunctions.CosineSeries(nbasis)
    Vout6 = zeros(length(wavetarg),1,nbasis)
    for chipind = 1:1
        mski = msklst[chipind]
        cntmski = length(mski)
        svec0 = range(0,1,length=cntmski)
        for (ind, sval) in enumerate(svec0)
            Vout6[mski[ind],chipind,:] .= cbasis(sval)
        end
    end

    Vreshape = reshape(hcat(Vout,Vout1[:,2],Vout2[:,2:end],Vout3[:,2],Vout4[:,2],Vout5[:,2],Vout6[:,1,2:4],Vout7[:,2]),length(wavetarg),:);

    return Vreshape, msknall
end

function ret_qlines(dat, waveobs;w = 500,rscale = 2)
    p = (0.25, 0.5, 0.75)
    moving_quartiles = movsort(dat, w, p)
    moving_q25 = circshift(map( x->x[1], moving_quartiles),-w÷2)
    moving_q50 = circshift(map( x->x[2], moving_quartiles),-w÷2)
    moving_q75 = circshift(map( x->x[3], moving_quartiles),-w÷2);

    qhigh = moving_q50 .+ rscale .*(moving_q75.-moving_q50);
    qlow = moving_q50 .+ rscale .*(moving_q25.-moving_q50)
    tmsk = (qlow .< dat .<qhigh)
    return tmsk
end

function find_nearest_nz(x,ind;totlen=300)
    if x[ind] !=0
        return ind
    else
        for di = 1:299
            nind = ind + di
            if (0 < nind <= totlen)
                if (x[nind] !=0)
                    return nind
                end
            end
            nind = ind - di
            if (0 < nind <= totlen)
                if (x[nind] !=0)
                    return nind
                end
            end
        end
    end
    return NaN
end

function expand_msk(msk;rad=1)
    lmsk = length(msk)
    msknew = ones(Bool,lmsk)
    for i=1:length(msk)
        lindx = maximum((1,i-rad))
        rindx = minimum((i+rad,lmsk))
        mskslice = view(msk,lindx:rindx)
        if any(.!mskslice)
            msknew[i] = false
        end
    end
    return msknew
end

function lorentzian(x,amp,x0,Γ)
    Γr = Γ/2
    return amp*(Γr/π)/((x-x0)^2 + Γr^2)
end

function redden_mult(wave,Av,Rv)
    return 10 .^(-0.4 .*Av*CCM89(Rv=Rv).(wave))
end

function blackbody(T, λ; hplanck_cgs = 6.62607015e-27, c_cgs = 2.99792458e10, kboltz_cgs = 1.380649e-16 )
    h = hplanck_cgs # erg*s
    c = c_cgs #cm/s
    k = kboltz_cgs # erg/K

    2*h*c^2/λ^5 * 1/(exp(h*c/λ/k/T) - 1)
end

function find_zero_ranges(vector)
    zero_ranges = Vector{Tuple{Int, Int}}()  # Stores start and end indices of zero ranges
    in_zero_range = false
    start_index = 0

    for (index, value) in enumerate(vector)
        if value == 0
            if !in_zero_range
                in_zero_range = true
                start_index = index
            end
        else
            if in_zero_range
                push!(zero_ranges, (start_index, index - 1))
                in_zero_range = false
            end
        end
    end

    # Check for an open zero range at the end of the vector
    if in_zero_range
        push!(zero_ranges, (start_index, length(vector)))
    end

    return zero_ranges
end

function gauss1d_ew(ew,x0,sigma,pixcoord)
    return ew .* exp.(-0.5 .*((pixcoord.-x0)./sigma).^2)./(sqrt(2*pi*(sigma^2)))
end

function gauss1d(amp,x0,sigma,pixcoord)
    return amp .* exp.(-0.5 .*((pixcoord.-x0)./sigma).^2)./(sqrt(2*pi))
end

function gauss1d_deriv(amp,x0,sigma,pixcoord)
    return amp .* (pixcoord.-x0) ./(sigma^2) .* exp.(-0.5 .*((pixcoord.-x0)./sigma).^2)./(sqrt(2*pi))
end

function gauss1d_2deriv(amp,x0,sigma,pixcoord)
    return - amp .* ((pixcoord.-x0).^2 ./(sigma^4) .- 1 ./(sigma^2)) .* exp.(-0.5 .*((pixcoord.-x0)./sigma).^2)./(sqrt(2*pi))
end

function SLURM_prune_workers_per_node(num_workers_per_node;reduce_head_by=1)
    getinfo_worker(workerid::Int) = @getfrom workerid myid(), ThreadPinning.sched_getcpu(), gethostname()
    dfo = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]);
    head_host = gethostname()
    activateout = @capture_out begin
        idlst = getinfo_worker.(workers()); df = DataFrame(workerid=Int[],physcpu=Int[],hostname=String[]); push!(df,idlst...)
        gdf = groupby(df,:hostname)
        for sgdf in gdf
            workers2keep = num_workers_per_node
            if head_host .== sgdf.hostname[1]
                workers2keep -= reduce_head_by
            end
            rmprocs(sgdf.workerid[(workers2keep+1):end])
        end
        idlst = getinfo_worker.(workers()); push!(dfo,idlst...)
    end
    # Helpful Worker Info Printing
    gdf = groupby(dfo,:hostname); dfc = combine(gdf, nrow, :workerid => minimum, :workerid => maximum, :physcpu => minimum, :physcpu => maximum)
    println("$(gethostname()) running Main")
    for row in Tables.namedtupleiterator(dfc)
        println("$(row.hostname) running $(row.nrow) workers: $(row.workerid_minimum)->$(row.workerid_maximum)")
    end
    flush(stdout)
end