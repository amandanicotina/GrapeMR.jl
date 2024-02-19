using Test

@safetestset "Bx: True Control Field Gradient vs Finite Difference" begin
    using GrapeMR

    # Defining parameters
    N   = 1000;
    α   = π/2;
    t_c = 1.0;
    B1  = α/(γ_¹H*t_c);
    B1x, B1y = B1*ones(1, N), B1*ones(1, N);   
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    field_tests = ControlFields(γ_¹H*B1x, γ_¹H*B1y, 0.0, 0.0, t_c, B0, ΔB0)
    spin_test  = Spins([0.0; 0.0; 1.0], 0.5, 0.3, 0.0, "max")
    mag_test   = forward_propagation(field_tests, spin_test)
    iso_test   = Magnetization((mag_test,), (spin_test,))
    cost       = "Target One Spin"
    cost_func  = cost_functions[cost]

    # Finite difference
    Δcf   = 1e-6;
    fd_cf = finite_difference_field(cost_func, field_tests, spin_test, Δcf)

    # True Gradient
    true_grad = -get_gradient(field_tests, iso_test, Ix, cost)
    
    plot(fd_cf')
    plot!(true_grad')

    @test round.(fd_cf, digits=2) == round.(true_grad, digits=2)

end

@safetestset "By: True Control Field Gradient vs Finite Difference" begin
    using GrapeMR

    # Defining parameters
    N   = 1000;
    α   = π/2;
    t_c = 1.0;
    B1  = α/(γ_¹H*t_c);
    B1x, B1y = B1*ones(1, N), B1*ones(1, N);   
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    field_test = ControlFields(γ_¹H*B1x, γ_¹H*B1y, 0.0, 0.0, t_c, B0, ΔB0)
    spin_test  = Spins([0.0; 0.0; 1.0], 0.5, 0.3, 0.0, "max")
    mag_test   = forward_propagation(field_test, spin_test)
    iso_test   = Magnetization((mag_test,), (spin_test,))
    cost       = "Target One Spin"
    cost_func  = cost_functions[cost]

    # Finite difference
    Δcf   = 1e-6;
    fd_cf = finite_difference_field(cost_func, field_test, spin_test, Δcf)

    # True Gradient
    true_grad = -get_gradient(field_test, iso_test, Iy, cost)
    
    plot(fd_cf')
    plot!(true_grad')

    @test round.(fd_cf, digits=2) == round.(true_grad, digits=2)

end

