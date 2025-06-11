using Test
using GrapeMR

@testset "process_spin!" begin
    using GrapeMR
    using Test

    # Mock inputs
    n = 10
    B1ref = ones(Float64, 1, n)
    control_field = spline_RF(n, 1.0, B1ref)

    spin = Spin(
        [0.0, 0.0, 1.0], 
        [1000.0], [100.0], 
        [0.0], 0.0, [1.0], ["test"]
    )

    grape_params = GrapeParams(
        n,
        (iso) -> (sum(abs2, iso.magnetization.m), [0.0, 0.0, 0.0, 0.0]),
        Dict("B1x" => true, "B1y" => true)
    )

    params = Parameters(grape_params, OptimizationParams(1.0, 1, 1))
    output = GrapeOutput(Isochromat[], deepcopy(control_field), zeros(Float64, 1), params)

    ∇x = zeros(Float64, 1, n)
    ∇y = zeros(Float64, 1, n)
    mag = zeros(Float64, 4, n + 1)
    adj = zeros(Float64, 4, n + 1)

    process_spin!(∇x, ∇y, mag, adj, spin, control_field, grape_params, 1, output)

    @test output.cost_values[1] > 0
    @test sum(abs, ∇x) > 0
    @test sum(abs, ∇y) > 0
end

@testset "grape!" begin
    using GrapeMR
    using Test

    n = 10
    control_field = spline_RF(n, 1.0, ones(Float64, 1, n))

    spin = Spin(
        [0.0, 0.0, 1.0],
        [1000.0], [100.0],
        [0.0], 0.0, [1.0], ["test"]
    )

    grape_params = GrapeParams(
        n,
        (iso) -> (sum(abs2, iso.magnetization.m), [0.0, 0.0, 0.0, 0.0]),
        Dict("B1x" => true, "B1y" => true)
    )

    opt_params = OptimizationParams(1.0, 1, 3)
    params = Parameters(grape_params, opt_params)

    output = GrapeOutput(Isochromat[], deepcopy(control_field), zeros(Float64, opt_params.max_iter), params)

    grape!(output, params, control_field, [spin])

    @test length(output.cost_values) == opt_params.max_iter
    @test all(output.cost_values .>= 0)
    @test length(output.isochromats) == 1
end

@testset "grape (wrapper)" begin
    using GrapeMR
    using Test

    n = 5
    B1ref = ones(Float64, 1, n)
    control_field = spline_RF(n, 1.0, B1ref)

    spin = Spin(
        [0.0, 0.0, 1.0],
        [1000.0], [100.0],
        [0.0], 0.0, [1.0], ["test"]
    )

    grape_params = GrapeParams(
        n,
        (iso) -> (sum(abs2, iso.magnetization.m), [0.0, 0.0, 0.0, 0.0]),
        Dict("B1x" => true, "B1y" => true)
    )

    opt_params = OptimizationParams(1.0, 1, 2)
    params = Parameters(grape_params, opt_params)

    out = grape(params, control_field, [spin])

    @test out isa GrapeOutput
    @test out.params == params
    @test out.control_field !== control_field # must be a deepcopy
end
