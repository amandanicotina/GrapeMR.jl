
using Test
using GrapeMR

@testset "Bloch Methods - Normalized ControlField" begin
    t_c = 0.01
    B1ref = 4.0
    cf = generate_control_field(:hard; t_c=t_c, B1ref=B1ref)  # Normalized by default

    spins = generate_spins(
        [0.0, 0.0, 1.0],
        [1.0], [0.5],
        [-10.0, 0.0, 10.0],
        [0.8, 1.0, 1.2],
        ["z", "z", "z"],
        ["test", "test", "test"]        
    )

    for spin in spins
        N = size(cf.B1x, 2)
        M = zeros(4, N + 1)
        forward_propagation!(M, cf, spin)

        @test size(M) == (4, N + 1)
        @test all(-1 .≤ M[2:4, :] .≤ 1)

        χ = zeros(4, N + 1)
        grad = [0.0, 0.0, 0.0, 1.0]
        iso = Isochromat(GrapeMR.Magnetization(M), spin)

        backward_propagation!(χ, cf, iso, grad)

        @test size(χ) == (4, N + 1)
        @test all(isfinite, χ)
    end
end

@testset "Bloch Methods - SI Units ControlField" begin
    N = 100
    t_c = 0.01
    B1ref = 4.0
    cf_si = generate_control_field(:hard; t_c=t_c, B1ref=B1ref, normalize=false)

    spins = generate_spins(
        [0.0, 0.0, 1.0],
        [1.0], [0.5],
        [-10.0, 0.0, 10.0],
        [0.8, 1.0, 1.2],
        ["z", "z", "z"],
        ["test", "test", "test"]        
    )

    for spin in spins
        N = size(cf_si.B1x, 2)
        M_si = zeros(4, N + 1)
        forward_propagation!(M_si, cf_si, spin)

        @test size(M_si) == (4, N + 1)
        @test all(-1 .≤ M_si[2:4, :] .≤ 1)

        χ_si = zeros(4, N + 1)
        grad_si = [0.0, 0.0, 0.0, 1.0]
        iso_si = Isochromat(GrapeMR.Magnetization(M_si), spin)

        backward_propagation!(χ_si, cf_si, iso_si, grad_si)

        @test size(χ_si) == (4, N + 1)
        @test all(isfinite, χ_si)
    end
end
