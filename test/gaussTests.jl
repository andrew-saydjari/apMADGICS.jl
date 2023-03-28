@testset "gridSearch.jl" begin
    
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
    x_comp_lst = deblend_components_all_asym_tot(Ctotinv, Xd_obs, 
        (A, Vgauss_r),
        (A, Vgauss_r),
    );
    
    dat2plot1 = nanify(x_data[simplemsk],simplemsk)
    dat2plot2 = nanify(x_comp_lst_asym[2][simplemsk],simplemsk)
    dat2plot3 = nanify(x_comp_lst[1],simplemsk)
    
    @test maximum(abs.(filter(.!isnan,(dat2plot2.+dat2plot3).-dat2plot1))) < 1e-12
    @test abs(sum(filter(.!isnan,(dat2plot2.+dat2plot3).-dat2plot1))) < 1e-10
    @test abs(sum(x_comp_lst_asym[2].-x_data_clean)) < 1e-2
    
end
