@testset "gridSearch.jl" begin
    ## 1D
    tstfun(x) = (x-1.5)^2+3
    tstrng = 0:5

    out = sampler_1d_dense(tstfun,tstrng)
    @test out[1] == (1.5, 3.0, 1, 3.25, 2, 0)

    out = sampler_1d_dense_var(tstfun,tstrng)
    @test out[1] == (1.5, 3.0, 1, 3.25, 2, 0, 0.5)

    tstrng1 = 0:2:6
    tstrng2 = -1:1
    tsttup = (tstrng1,tstrng2)

    out = sampler_1d_hierarchy_var(tstfun,tsttup)
    @test out[1] == (1.5, 3.0, 3//2, 3.0, 2, 0, 0.5)
    
    ## 2D
    function tstfun2d(svals)
        x, y = svals
        return (x-1.5)^2 + (y-2.5)^2 + 3
    end
    tstrng_2d = ((0:5),(0:5));
    
    out = sampler_2d_dense(tstfun2d,tstrng_2d)
    @test out[1] == ((1.5, 2.5), 3.0, (1, 2), 3.5, CartesianIndex(2, 3), 0)
    
    out = sampler_2d_dense_var(tstfun2d,tstrng_2d)
    @test out[1] == ((1.5, 2.5), 3.0, (1, 2), 3.5, CartesianIndex(2, 3), 0, 0.5, 0.5, -0.0, 0.5, 0.5)
    
    tstrng_2d_1 = ((0:5),(0:5));
    tstrng_2d_2 = ((-2:2),(-2:2));
    tsttup_2d = (tstrng_2d_1,tstrng_2d_2);
    
    out = sampler_2d_hierarchy_var(tstfun2d,tsttup_2d)
    @test out[1] == ((1.5, 2.5), 3.0, (3//2, 5//2), 3.0, CartesianIndex(3, 3), 0, 0.5, 0.5, -0.0, 0.5, 0.5)
    
end
