using Test
using GrapeMR
# using Plots  # Uncomment if you want to visualize

@testset "Gradient Check: Adjoint vs Finite Differences (B1x & B1y)" begin
    # Parameters
    t_c = 0.1
    B1ref = 4.0
    δ = 1e-5

    # Setup control field and spin
    cf = generate_control_field(:spline; t_c=t_c, B1ref=B1ref)
    N = length(cf.B1x)

    m_init = [0.0, 0.0, 1.0]
    spin = Spin(m_init, 1.0, 0.2, 0.0, 1.0, "origin", "test", 1)

    # Forward propagation
    M = zeros(4, N + 1)
    forward_propagation!(M, cf, spin)
    iso = Isochromat(GrapeMR.Magnetization(M), spin)
    val, cost_grad = euclidean_norm(iso)

    # Backward propagation
    χ = zeros(4, N + 1)
    backward_propagation!(χ, cf, iso, cost_grad)

    # Adjoint gradients
    grad_x = zeros(1, N)
    grad_y = zeros(1, N)
    gradient!(grad_x, χ, M, Ix)
    gradient!(grad_y, χ, M, Iy)

    # Finite differences
    grad_fd_x = zeros(1, N)
    grad_fd_y = zeros(1, N)

    for i in 1:N
        cf_p = deepcopy(cf); cf_p.B1x[i] += δ
        M_p = zeros(4, N + 1)
        forward_propagation!(M_p, cf_p, spin)
        iso_p = Isochromat(GrapeMR.Magnetization(M_p), spin)
        val_p, _ = euclidean_norm(iso_p)
    
        cf_m = deepcopy(cf); cf_m.B1x[i] -= δ
        M_m = zeros(4, N + 1)
        forward_propagation!(M_m, cf_m, spin)
        iso_m = Isochromat(GrapeMR.Magnetization(M_m), spin)
        val_m, _ = euclidean_norm(iso_m)
    
        grad_fd_x[1, i] = (val_p - val_m) / (2δ)
    end    

    for i in 1:N
        cf_p = deepcopy(cf); cf_p.B1y[i] += δ
        M_p = zeros(4, N + 1)
        forward_propagation!(M_p, cf_p, spin)
        iso_p = Isochromat(GrapeMR.Magnetization(M_p), spin)
        val_p, _ = euclidean_norm(iso_p)
    
        cf_m = deepcopy(cf); cf_m.B1y[i] -= δ
        M_m = zeros(4, N + 1)
        forward_propagation!(M_m, cf_m, spin)
        iso_m = Isochromat(GrapeMR.Magnetization(M_m), spin)
        val_m, _ = euclidean_norm(iso_m)
    
        grad_fd_y[1, i] = (val_p - val_m) / (2δ)
    end
    

    # Diagnostics
    max_Δx = maximum(abs.(grad_x .- grad_fd_x))
    rel_err_x = max_Δx / maximum(abs.(grad_x))

    max_Δy = maximum(abs.(grad_y .- grad_fd_y))
    rel_err_y = max_Δy / maximum(abs.(grad_y))

    println("Max Δ B1x: ", max_Δx, " | Relative Error: ", rel_err_x)
    println("Max Δ B1y: ", max_Δy, " | Relative Error: ", rel_err_y)

    # Optional plotting (uncomment if needed)
    # plot(grad_x[:], label="Adjoint ∇B₁x")
    # plot!(grad_fd_x[:], label="FD ∇B₁x", linestyle=:dash)
    # title!("Gradient Check: B1x")

    # plot(grad_y[:], label="Adjoint ∇B₁y")
    # plot!(grad_fd_y[:], label="FD ∇B₁y", linestyle=:dash)
    # title!("Gradient Check: B1y")

    # Assertions
    @test isapprox(grad_x, grad_fd_x; rtol=5e-2, atol=1e-4)
    @test isapprox(grad_y, grad_fd_y; rtol=5e-2, atol=1e-4)

end
