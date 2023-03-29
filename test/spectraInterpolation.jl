@testset "spectraInterpolation.jl" begin
  
  wavetarg = 10 .^range((4.179-125*6.0e-6),step=6.0e-6,length=8575+125);
  testset = [15000.0,15273.0,15263.0,17000.0]
  function checkfindmin(x)
    return findmin(abs.(wavetarg.-x))[2]
  end
  targind = find_yinx(wavetarg,testset)
  refind = checkfindmin.(testset)
  
  @test targind == refind 
end
