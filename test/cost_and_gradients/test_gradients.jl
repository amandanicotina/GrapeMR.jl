using Test
using GrapeMR

@testset "Gradient calculation and update" begin
    # Setup
    t_c = 0.01
    B1ref = 4.0
    cf = generate_control_field(:hard; t_c=t_c, B1ref=B1ref)  # Normalized by default

    # Spin setup
    m_init = [0.0, 0.0, 1.0]
    spin = Spin(m_init, 1.0, 0.2, 0.0, 1.0, "max", "test", 1)

    # Forward propagation
    N = size(cf.B1x, 2)
    M = zeros(4, N + 1)
    forward_propagation!(M, cf, spin)

    # Isochromat
    iso = Isochromat(GrapeMR.Magnetization(M), spin)

    # Cost function and gradient
    val, cost_grad = saturation_contrast(iso)
    @test val isa Number
    @test length(cost_grad) == 4

    # Backward propagation
    χ = zeros(4, N + 1)
    backward_propagation!(χ, cf, iso, cost_grad)
    @test size(χ) == (4, N + 1)

    grad_x = zeros(1, N)
    grad_y = zeros(1, N)

    gradient!(grad_x, χ, M, Ix)
    gradient!(grad_y, χ, M, Iy)

    @test size(grad_x) == (1, N)
    @test size(grad_y) == (1, N)

    # Update
    ϵ = 0.1
    u1x, u1y = update!(cf, (grad_x, grad_y), ϵ)
    @test length(u1x) == N
    @test length(u1y) == N
end









