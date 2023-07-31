@testset "utils.jl" begin

    dat = 2.1*ones(10,10)
    dat[1,:].=NaN

    @test nanmean(dat[:,1])==2.1
    @test nanmean(dat,1)==2.1*ones(1,10)

    @test nansum(dat[:,1],1)[1] == 2.1*9
    @test nansum(dat,1) == 2.1*9*ones(1,10)

    @test nanmedian(dat[:,1])==2.1
    @test nanmedian(dat,1)==2.1*ones(1,10)

    # Gaussian Posterior Test
    using Distributions
    x = 0.5
    x0 = 0
    s0 = 1
    # Calculate the expected value using Distributions.jl
    d = Normal(x0, s0)
    expected = pdf(d, x)
    @test isapprox(gaussian_post(x, x0, s0), expected)

    x = 3.5
    x0 = 1
    s0 = 2
    # Calculate the expected value using Distributions.jl
    d = Normal(x0, s0)
    expected = pdf(d, x)
    @test isapprox(gaussian_post(x, x0, s0), expected)

    # Test for positive input
    @test isapprox(sqrt_nan(9), 3.0, rtol=1e-6)
    @test isapprox(sqrt_nan(0), 0.0, rtol=1e-6)

    # Test for negative input
    @test isnan(sqrt_nan(-5))

    # Test for NaN input
    @test isnan(sqrt_nan(NaN))

    # Test for edge case of Inf input
    @test isinf(sqrt_nan(Inf))

    # Test for a square matrix (2x2)
    @test issquare([1 2; 3 4])

    # Test for a non-square matrix (2x3)
    @test !issquare([1 2 3; 4 5 6])

    # Covariance Ellipse Test
    Σ = [1.0 0.0; 0.0 1.0]
    μ = [0.0, 0.0]
    n_std = 1
    n_ellipse_vertices = 100
    # Calculate ellipse points using covellipse
    ellipse_points = covellipse(Σ, μ=μ, n_std=n_std, n_ellipse_vertices=n_ellipse_vertices)
    
    # Calculate ellipse points using the formula for an ellipse with identity covariance matrix
    θ = range(0, 2π; length=n_ellipse_vertices)
    x = μ[1] .+ n_std * cos.(θ)
    y = μ[2] .+ n_std * sin.(θ);
    @test isapprox(ellipse_points[1], x)
    @test isapprox(ellipse_points[2], y)

    Σ = [4.0 0.0; 0.0 4.0]
    μ = [0.0, 0.0]
    n_std = 1
    n_ellipse_vertices = 100
    # Calculate ellipse points using covellipse
    ellipse_points = covellipse(Σ, μ=μ, n_std=n_std, n_ellipse_vertices=n_ellipse_vertices)
    
    # Calculate ellipse points using the formula for an ellipse with identity covariance matrix
    θ = range(0, 2π; length=n_ellipse_vertices)
    x = μ[1] .+ n_std * 2.0 * cos.(θ)
    y = μ[2] .+ n_std * 2.0 *sin.(θ);
    @test isapprox(ellipse_points[1], x)
    @test isapprox(ellipse_points[2], y)

    lvl1 = -70:1//2:70
    lvl2 = -8:2//10:8
    lvl3 = -3:1//10:3
    slvl_tuple = (lvl1,lvl2,lvl3)
    @test tuple1dprint(slvl_tuple) == println([281, 81, 61, 423])

    lvl1d = ((-150:4:150),(18//10:18//10))
    lvl2d = ((0:0), (-7//5:4//100:11//5))
    lvl3d = ((-18:2//10:18), (0:0))
    lvltuple = (lvl1d, lvl2d, lvl3d);
    @test tuple2dprint(lvltuple) == println([76, 91, 181, 348])

    x = [1.1, 2.5, 3.7, 4.2, 5.8]
    msk = [true, false, true, false, true]
    y = nanify(x[msk], msk)
    @test y[msk] == [1.1, 3.7, 5.8]
    @test all(isnan.(y[.!msk]))

    # cpu_lock test (rough)
    # @test_throws nothing slurm_cpu_lock()
    using ThreadPinning, DataFrames, Distributed
    slurm_cpu_lock()

    z = v2z(10)
    v = z2v(z)
    @test v ≈ 10

    v = pix2v(1)
    @test v ≈ 4.141785868417105
    z = v2z(v)
    p = z2pix(z)
    @test p ≈ 1

    @test prop_p2z(1) ≈ 2.302616904602455
    @test prop_z2v(1)==1.6
end