using Test
using GrapeMR

@testset "Cost Function Behavior: euclidean_norm, spin_target, saturation_contrast" begin
    # Common setup
    t_c = 0.01
    B1ref = 4.0
    cf = generate_control_field(:hard; t_c=t_c, B1ref=B1ref)  # Simple pulse

    m_init = [0.0, 0.0, 1.0]
    spin = Spin(m_init, 1.0, 0.5, 0.0, 1.0, "max", "test", 1)

    # Forward propagation
    N = size(cf.B1x, 2)
    M = zeros(4, N + 1)
    forward_propagation!(M, cf, spin)
    iso = Isochromat(GrapeMR.Magnetization(M), spin)

    ## euclidean_norm
    val_euclid, grad_euclid = euclidean_norm(iso)
    @test val_euclid isa Number
    @test length(grad_euclid) == 4
    @test grad_euclid[1] == 0.0  # gradient w.r.t time component should be zero

    ## spin_target — target = m_final
    m_final = M[2:4, end]  # use final magnetization as perfect target
    val_target, grad_target = spin_target(iso; target=collect(m_final))
    @test isapprox(val_target, 0.0; atol=1e-5)
    @test grad_target[1] == 0.0
    @test all(abs.(grad_target[2:4]) .< 1e-8)

    ## saturation_contrast — target = "max"
    spin_max = Spin(m_init, 1.0, 0.5, 0.0, 1.0, "max", "max_spin", 1)
    iso_max = Isochromat(GrapeMR.Magnetization(M), spin_max)
    val_sat_max, grad_sat_max = saturation_contrast(iso_max)
    @test val_sat_max isa Number
    @test length(grad_sat_max) == 4
    @test grad_sat_max[1] == 0.0

    ## saturation_contrast — target = "min"
    spin_min = Spin(m_init, 1.0, 0.5, 0.0, 1.0, "min", "min_spin", 1)
    iso_min = Isochromat(GrapeMR.Magnetization(M), spin_min)
    val_sat_min, grad_sat_min = saturation_contrast(iso_min)
    @test val_sat_min isa Number
    @test length(grad_sat_min) == 4
    @test grad_sat_min[1] == 0.0
end
