@testset "gaussTests.jl" begin
    
    ### Gaussian peak gap reconstruction
    function gauss1d(amp,x0,sigma,pixcoord)
        return amp .* exp.(-0.5 .*((pixcoord.-x0)./sigma).^2)./(sqrt(2*pi))
    end

    function gauss1d_deriv(amp,x0,sigma,pixcoord)
        return - amp .* (pixcoord.-x0) ./(sigma^2) .* exp.(-0.5 .*((pixcoord.-x0)./sigma).^2)./(sqrt(2*pi))
    end

    wavemodel = 15220:0.1:15320
    minw, maxw = extrema(wavemodel)
    
    simplemsk = ones(Bool,length(wavemodel))
    simplemsk[560:620] .= false

    Vgauss = zeros(length(wavemodel),2)
    cwave_cent = 15273
    amp = -0.1
    damp = 0.01
    offset = 0
    sigma = 1.9
    x0 = cwave_cent*10^(offset*6e-6)
    Vgauss[:,1] .= gauss1d(amp,x0,sigma,wavemodel)
    Vgauss[:,2] .= gauss1d_deriv(damp*amp,x0,sigma,wavemodel);

    amp = -0.1
    noiseFac = 1e-3
    rng = MersenneTwister(2025)
    var_vec = (noiseFac.*range(1,2,length(wavemodel))).^2
    x_data_clean = gauss1d(amp,x0,sigma,wavemodel)
    x_noise_vec = sqrt.(var_vec).*randn(rng,length(wavemodel));
    x_data = x_data_clean .+ x_noise_vec;
    
    ## Select data for use (might want to handle mean more generally)
    Xd_obs = x_data[simplemsk]; 
    wave_obs = wavemodel[simplemsk]

    ## Set up residuals prior
    A = Diagonal(var_vec[simplemsk]);
    Ainv = Diagonal(1 ./var_vec[simplemsk]);

    ## Set up priors
    Vgauss_c = Vgauss
    Vgauss_r = Vgauss_c[simplemsk,:]

    ## Solve RV of Star
    # compute stellar continuum to modify stellar line prior
    Ctotinv = LowRankMultMat([Ainv,Vgauss_r],wood_precomp_mult,wood_fxn_mult);
    x_comp_lst_asym = deblend_components_all_asym_tot(Ctotinv, Xd_obs, 
        (A, Vgauss_r),
        (A, Vgauss_c),
    )
    x_comp_lst = deblend_components_all_tot(Ctotinv, Xd_obs, 
        (A, Vgauss_r),
    );
    
    dat2plot1 = nanify(x_data[simplemsk],simplemsk)
    dat2plot2 = nanify(x_comp_lst_asym[2][simplemsk],simplemsk)
    dat2plot3 = nanify(x_comp_lst[1],simplemsk)
    
    @test maximum(abs.(filter(.!isnan,(dat2plot2.+dat2plot3).-dat2plot1))) < 1e-12
    @test abs(sum(filter(.!isnan,(dat2plot2.+dat2plot3).-dat2plot1))) < 1e-10
    @test abs(sum(x_comp_lst_asym[2].-x_data_clean)) < 1e-2
    
    ## Gaussian interpolation linear->log with missing pixels/chip gaps
    waveobs = 15220:0.1:15320
    wavemodel = 10 .^(log10(15220):3e-6:log10(15320))
    minw, maxw = extrema(wavemodel)
    amp = -0.1
    rng = MersenneTwister(2035)
    function std_fxn(x;noiseFac=5e-3)
        return noiseFac*(1+0.001*(x-15220)/(15320-15220))
    end
    var_vec = (std_fxn.(waveobs)).^2
    x_data_clean = gauss1d(amp,x0,sigma,waveobs)
    x_noise_vec = sqrt.(var_vec).*randn(rng,length(waveobs));
    x_data = x_data_clean .+ x_noise_vec;

    x_targ = gauss1d(amp,x0,sigma,wavemodel)

    simplemsk = ones(Bool,length(waveobs))
    simplemsk[560:620] .= false
    simplemsk[500:501] .= false;
    
    simplemsk[500:501] .= false;

    fullBit = zeros(Int,length(waveobs));
    fullBit[1:559] .+= 2^1 # theoretical chip1
    fullBit[621:end] .+= 2^2 # theoretical chip2
    fullBit[500:501] .+= 2^4 # bad pixels

    Xd_obs = x_data[simplemsk]; 
    wave_obs = waveobs[simplemsk]
    obsBit = fullBit[simplemsk]
    pixindx = (1:length(waveobs))[simplemsk]

    Rinv = generateInterpMatrix_sparse_inv(wave_obs,obsBit,wavemodel,pixindx)
    normvec = dropdims(sum(Rinv,dims=2),dims=2)
    propmsk = (normvec.!=0)
    
    @test maximum((Rinv*x_data_clean[simplemsk].-x_targ)[propmsk]) .< 3e-5
    @test sum((Rinv*x_data_clean[simplemsk].-x_targ)[propmsk]) .< 5e-5
    
    waveobs = (sqrt(15220):0.0007:sqrt(15320)).^2
    wavemodel = 10 .^(log10(15220):3e-6:log10(15320))
    minw, maxw = extrema(wavemodel)
    amp = -0.1
    rng = MersenneTwister(2035)
    function std_fxn(x;noiseFac=5e-3)
        return noiseFac*(1+0.001*(x-15220)/(15320-15220))
    end
    var_vec = (std_fxn.(waveobs)).^2
    x_data_clean = gauss1d(amp,x0,sigma,waveobs)
    x_noise_vec = sqrt.(var_vec).*randn(rng,length(waveobs));
    x_data = x_data_clean .+ x_noise_vec;

    x_targ = gauss1d(amp,x0,sigma,wavemodel)

    simplemsk = ones(Bool,length(waveobs))
    simplemsk[320:330] .= false
    simplemsk[305:305] .= false;

    fullBit = zeros(Int,length(waveobs));
    fullBit[1:319] .+= 2^1 # theoretical chip1
    fullBit[331:end] .+= 2^2 # theoretical chip2
    fullBit[305:305] .+= 2^4 # bad pixels

    Xd_obs = x_data[simplemsk]; 
    wave_obs = waveobs[simplemsk]
    obsBit = fullBit[simplemsk]
    pixindx = (1:length(waveobs))[simplemsk]

    Rinv = generateInterpMatrix_sparse_inv(wave_obs,obsBit,wavemodel,pixindx)
    normvec = dropdims(sum(Rinv,dims=2),dims=2)
    propmsk = (normvec.!=0);

    Xd_obs = (Rinv*x_data[simplemsk])[propmsk]; 
    wave_obs = wavemodel[propmsk]
    var_obs = ((Rinv.^2)*var_vec[simplemsk])

    Vgauss = zeros(length(wavemodel),2)
    cwave_cent = 15273
    amp = -0.1
    damp = 0.01
    offset = 0
    sigma = 1.9
    x0 = cwave_cent*10^(offset*6e-6)
    Vgauss[:,1] .= gauss1d(amp,x0,sigma,wavemodel)
    Vgauss[:,2] .= gauss1d_deriv(damp*amp,x0,sigma,wavemodel);

    ## Set up residuals prior
    A = Diagonal(var_obs[propmsk]);
    Ainv = Diagonal(1 ./var_obs[propmsk]);

    ## Set up priors
    Vgauss_c = Vgauss
    Vgauss_r = Vgauss_c[propmsk,:]

    ## Solve RV of Star
    # compute stellar continuum to modify stellar line prior
    Ctotinv = LowRankMultMat([Ainv,Vgauss_r],wood_precomp_mult,wood_fxn_mult);
    x_comp_lst_asym = deblend_components_all_asym_tot(Ctotinv, Xd_obs, 
        (A, Vgauss_r),
        (A, Vgauss_c),
    )
    x_comp_lst = deblend_components_all_tot(Ctotinv, Xd_obs, 
        (A, Vgauss_r),
    );

    dat2plot1 = nanify(Xd_obs,propmsk)
    dat2plot2 = nanify(x_comp_lst_asym[2][propmsk],propmsk)
    dat2plot3 = nanify(x_comp_lst[1],propmsk);
    
    @test maximum((Rinv*x_data_clean[simplemsk].-x_targ)[propmsk]) .< 5e-5
    @test sum((Rinv*x_data_clean[simplemsk].-x_targ)[propmsk]) .< 3e-4
    
    @test maximum(abs.(filter(.!isnan,(dat2plot2.+dat2plot3).-dat2plot1))) < 1e-12
    @test abs(sum(filter(.!isnan,(dat2plot2.+dat2plot3).-dat2plot1))) < 1e-10
    @test abs(sum(x_comp_lst_asym[2].-x_targ)) < 1e-2 
end
