@testset "spectraInterpolation.jl" begin
  
  wavetarg = 10 .^range(start=(4.179-125*6.0e-6),step=6.0e-6,length=8575+125);
  testset = [15000,15273,15263,17000]
  function checkfindmin(x)
    return findmin(abs.(wavetarg.-x))[2]
  end
  
  @test find_yinx(wavetarg,testset) == checkfindmin.(testset)
end
