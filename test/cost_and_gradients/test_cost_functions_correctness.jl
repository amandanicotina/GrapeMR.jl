using Test
using GrapeMR

@testset "Gradient Check: euclidean_norm" begin
    t_c = 0.1
    B1ref = 4.0
    ε = 1e-6

    cf = generate_control_field(:hard; t_c=t_c, B1ref=B1ref)
    spin = Spin([0.0, 0.0, 1.0], 1.0, 0.5, 0.0, 1.0, "norm_test", "euclidean", 1)
    
    N = size(cf.B1x, 2)
    M = zeros(4, N + 1)
    forward_propagation!(M, cf, spin)
    iso = Isochromat(GrapeMR.Magnetization(M), spin)

    val_analytical, grad_analytical = euclidean_norm(iso)

    grad_fd = zeros(length(grad_analytical))
    for i in 2:4
        m_dyn = copy(iso.magnetization.dynamics)

        m_dyn[i, end] += ε
        iso_plus = Isochromat(GrapeMR.Magnetization(m_dyn), spin)
        val_plus, _ = euclidean_norm(iso_plus)

        m_dyn[i, end] -= 2ε
        iso_minus = Isochromat(GrapeMR.Magnetization(m_dyn), spin)
        val_minus, _ = euclidean_norm(iso_minus)

        grad_fd[i] = (val_plus - val_minus) / (2ε)
    end

    @test isapprox(grad_fd[2:4], grad_analytical[2:4]; atol=1e-4)
end

@testset "Gradient Check: spin_target" begin
    t_c = 0.01
    B1ref = 4.0
    ε = 1e-6

    cf = generate_control_field(:hard; t_c=t_c, B1ref=B1ref)
    target_vec = [0.0, 1.0, 0.0]

    spin = Spin([0.0, 0.0, 1.0], 1.0, 0.5, 0.0, 1.0, "target_test", "spin_target", 1)

    N = size(cf.B1x, 2)
    M = zeros(4, N + 1)
    forward_propagation!(M, cf, spin)
    iso = Isochromat(GrapeMR.Magnetization(M), spin)

    val_analytical, grad_analytical = spin_target(iso; target=target_vec)

    grad_fd = zeros(length(grad_analytical))
    for i in 2:4
        m_dyn = copy(iso.magnetization.dynamics)

        m_dyn[i, end] += ε
        iso_plus = Isochromat(GrapeMR.Magnetization(m_dyn), spin)
        val_plus, _ = spin_target(iso_plus; target=target_vec)

        m_dyn[i, end] -= 2ε
        iso_minus = Isochromat(GrapeMR.Magnetization(m_dyn), spin)
        val_minus, _ = spin_target(iso_minus; target=target_vec)

        grad_fd[i] = (val_plus - val_minus) / (2ε)
    end

    @test isapprox(grad_fd[2:4], grad_analytical[2:4]; atol=1e-4)
end


@testset "Gradient Check: saturation_contrast" begin
    t_c = 0.01
    B1ref = 4.0
    ε = 1e-6

    cf = generate_control_field(:hard; t_c=t_c, B1ref=B1ref)

    # test "max" target case
    spin = Spin([0.0, 0.0, 1.0], 1.0, 0.5, 0.0, 1.0, "max", "contrast_max", 1)
    
    N = size(cf.B1x, 2)
    M = zeros(4, N + 1)
    forward_propagation!(M, cf, spin)
    iso = Isochromat(GrapeMR.Magnetization(M), spin)

    val_analytical, grad_analytical = saturation_contrast(iso)

    grad_fd = zeros(length(grad_analytical))
    for i in 2:4
        m_dyn = copy(iso.magnetization.dynamics)

        m_dyn[i, end] += ε
        iso_plus = Isochromat(GrapeMR.Magnetization(m_dyn), spin)
        val_plus, _ = saturation_contrast(iso_plus)

        m_dyn[i, end] -= 2ε
        iso_minus = Isochromat(GrapeMR.Magnetization(m_dyn), spin)
        val_minus, _ = saturation_contrast(iso_minus)

        grad_fd[i] = (val_plus - val_minus) / (2ε)
    end

    @test isapprox(grad_fd[2:4], grad_analytical[2:4]; atol=1e-4)
end
