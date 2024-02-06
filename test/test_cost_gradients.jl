using Test

@safetestset "True Cost Function Gradient vs Finite Difference" begin
    using GrapeMR
    # Defining parameters
    N = 100;
    fields_test = InitialControlFields(N, rand(1, N), 1.0, 0.0, 1e-1, 0.0, 0.0)
    spin_test   = Spins([1.0; 0.0; 0.0], 0.5, 0.3, 0.0, "max")
    mag_test    = magnetization_ODE(fields_test, spin_test)
    iso_test    = Magnetization((mag_test,), (spin_test,))

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = cost_functions["Euclidean Norm"]

    # Finite difference
    ΔM = 1e-8;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_gradients["Grad Euclidean Norm"](iso_test)[2:end,:]

    @test round.(fd_M, digits=5) == round.(true_grad, digits=5)
end


@safetestset "True Control Field Gradient vs Finite Difference" begin
    using GrapeMR

    # Defining parameters
    N = 500;
    α   = π/2;
    t_c = 1.0;
    B1  = α/(γ_¹H*t_c);
    B1x_arr, B1y_arr = B1*ones(1, N), zeros(1,N);   

    fields_test = InitialControlFields(N, B1x_arr, 1.0, B1y_arr, 0.0, 1e-1, 0.0, 0.0)
    spin_test   = Spins([1.0; 0.0; 0.0], 0.5, 0.3, 0.0, "max")
    mag_test    = magnetization_ODE(fields_test, spin_test)
    iso_test    = Magnetization((mag_test,), (spin_test,))
    cost_func   = cost_functions["Euclidean Norm"]

    # Finite difference
    Δcf   = 1e-10;
    fd_cf = finite_difference_field(cost_func, fields_test, spin_test, Δcf)./γ_¹H

    # True Gradient
    true_grad = -gradient_controls(fields_test, spin_test, iso_test)
plot(fd_cf')
plot!(true_grad')

    @test round.(fd_cf, digits=5) .== round.(true_grad, digits=5)

end

