# tests/test_generate_control_field.jl
using Test
using GrapeMR
using Random

@testset "ControlField construction and normalization" begin
    N = 100
    t_c = 0.01  # 10 ms
    B1ref = 5.0

    B1x = B1ref .* ones(N)
    B1y = zeros(N)
    Bz = 0.5 .* ones(N)

    cf = GrapeMR.build_control_field(B1x, B1y, B1ref, t_c, Bz=Bz)

    @test cf.B1x isa Matrix
    @test size(cf.B1x) == (1, N)
    @test cf.B1_ref ≈ B1ref
    @test cf.t_control ≈ t_c

    cf_norm = GrapeMR.normalize_control_field(cf; t_unit=1e-3, B1_unit=B1ref)
    @test cf_norm.B1_ref ≈ 1.0
    @test cf_norm.t_control ≈ 10.0  # 10 ms in ms

    cf_back = GrapeMR.denormalize_control_field(cf_norm; t_unit=1e-3, B1_unit=B1ref)
    @test cf_back.B1_ref ≈ B1ref
    @test cf_back.t_control ≈ t_c
    @test all(abs.(cf_back.B1x .- cf.B1x) .< 1e-10)
end

@testset "Pulse generation: all types normalized" begin
    N = 100
    t_c = 0.01
    B1ref = 4.0

    for kind in (:spline, :hard, :sinc, :gaussian)
        cf = generate_control_field(kind; N=N, t_c=t_c, B1ref=B1ref)

        @test cf isa NormalizedControlField
        @test size(cf.B1x, 2) == N

        # Check magnitude constraint
        constrained_magnitude = sqrt.(cf.B1x.^2 .+ cf.B1y.^2)
        @test maximum(constrained_magnitude) ≤ 1.0
    end

    # bSSFP case
    cf_bssfp = generate_control_field(:bssfp; N=100, t_c=t_c, B1ref=B1ref, nTR=5, α=π/2, TR=0.002)
    @test cf_bssfp isa NormalizedControlField
    @test size(cf_bssfp.B1x, 2) == 100
end


@testset "Spline interpolation sanity check" begin
    rng = MersenneTwister(42)
    knots = range(0, 1, length=10)
    vals = rand(rng, 10)
    eval_pts = range(0, 1, length=100)
    spline_vals = GrapeMR.create_spline(knots, eval_pts, vals)

    @test length(spline_vals) == 100
    @test all(isfinite, spline_vals)
end
