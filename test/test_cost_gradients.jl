using Test

@safetestset "True Control Field Gradient vs Finite Difference" begin
    using GrapeMR

    # Defining parameters
    N = 500;
    α   = π/2;
    t_c = 1.0;
    B1  = α/(γ_¹H*t_c);
    B1x_arr, B1y_arr = B1*ones(1, N), zeros(1,N);   

    fields_test = InitialControlFields(N, γ_¹H*B1x_arr, 0.0, B1y_arr, 0.0, t_c, 0.0, 0.0)
    spin_test   = Spins([1.0; 0.0; 0.0], 0.5, 0.3, 0.0, "max")
    mag_test    = forward_propagation(fields_test, spin_test)
    iso_test    = Magnetization((mag_test,), (spin_test,))
    cost_func   = cost_functions["Target One Spin"]

    # Finite difference
    Δcf   = 1e-6;
    fd_cf = finite_difference_field(cost_func, fields_test, spin_test, Δcf)

    # True Gradient
    true_grad = -gradient_controls(fields_test, spin_test, iso_test)
    
    #plot(fd_cf')
    #plot!(true_grad')

    @test round.(fd_cf, digits=2) == round.(true_grad, digits=2)

end

