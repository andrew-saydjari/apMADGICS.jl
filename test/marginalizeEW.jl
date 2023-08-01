@testset "spectraInterpolation.jl" begin
  
  chi2lst = [1.0, 2.0, 3.0]
  fluxlst = [NaN, NaN, NaN]
  dfluxlst = [NaN, NaN, NaN]
  refchi2val = 1.0
  result = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
  @test length(result) == 2
  @test isnan(result[1])  # Check that the mean is NaN
  @test isnan(result[2])  # Check that the standard deviation is NaN

  chi2lst = [1.0, 1.0]
  fluxlst = [4.0, -4.0]
  dfluxlst = [1.0, 1.0]
  refchi2val = 1.0
  result = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
  @test abs(result[1]) < 1e-15

  chi2lst = [1.0]
  fluxlst = [5.0]
  dfluxlst = [2.0]
  refchi2val = 1.0
  result = marginalize_flux_err(chi2lst, fluxlst, dfluxlst, refchi2val)
  @test abs(result[1].-5)<1e-14
  @test abs(result[2].-2)<1e-4
  
end
