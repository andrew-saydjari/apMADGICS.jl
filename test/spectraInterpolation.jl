@testset "spectraInterpolation.jl" begin

  using StableRNGs
  
  wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125);
  testset = [15000.0,15273.0,15263.0,17000.0]
  function checkfindmin(x)
    return findmin(abs.(wavetarg.-x))[2]
  end
  targind = find_yinx(wavetarg,testset)
  refind = checkfindmin.(testset)
  
  @test targind == refind 

  # Testing that diagonal noise comes through interpolation correctly
  wave_line = sqrt(15000):(2//1000):sqrt(17000)
  wave_quad1 = wave_line.^2;
  wave_quad2 = wave_quad1 .+ StatsBase.median(diff(wave_quad1))/2;

  outvec = zeros(length(wavetarg))
  outvar = zeros(length(wavetarg))
  cntvec = zeros(Int,length(wavetarg));
  
  mu = 2
  sigin = 5
  exptot = 64
  rng = StableRNG(1486)
  fill!(outvec,0)
  fill!(outvar,0)
  fill!(cntvec,0)
  for i in 1:exptot
      waveobs = if i % 2 == 0 
          wave_quad1
      else
          wave_quad2
      end
      Xd_std_obs = sigin*ones(length(waveobs))
      Xd_obs = sigin*randn(rng,length(waveobs)).+mu;
      Rinv = generateInterpMatrix_sparse_inv(waveobs,ones(Int,length(waveobs)),wavetarg,1:length(waveobs))
      normvec = dropdims(sum(Rinv,dims=2),dims=2)
      msk_inter = (normvec.!=0)

      fullvec = Rinv*Xd_obs
      fullvec[.!msk_inter] .= 0
      varvec =  (Rinv.^2)*(Xd_std_obs.^2)
      varvec[.!msk_inter] .= 0

      outvec .+= fullvec
      outvar .+= varvec
      cntvec .+= msk_inter;
  end
  framecnts = maximum(cntvec) # a little shocked that I throw it away if it is bad in even one frame
  outvec./=framecnts
  outvar./=(framecnts^2);

  mu_est = nanzeromedian(outvec)
  z_est = nanzeroiqr((outvec.-mu_est)./sqrt.(outvar))

  @test abs(mu_est-mu) < 2e-3
  @test abs(z_est-1) < 1e-2

end
