using Test

@safetestset "True Cost Function Gradient vs Finite Difference" begin
    using GrapeMR
    # Defining parameters
    N = 500;
    α   = π/2;
    t_c = 0.1;
    B1  = α/(γ_¹H*t_c);
    B1x_arr, B1y_arr = B1*ones(1, N), zeros(1,N);   

    fields_test = InitialControlFields(N, γ_¹H*B1x_arr, 0.0, B1y_arr, 0.0, t_c, 0.0, 0.0)
    spin_test   = Spins([0.0; 1.0; 0.0], 0.5, 0.3, 0.0, "max")
    mag_test    = forward_propagation(fields_test, spin_test)
    iso_test    = Magnetization((mag_test,), (spin_test,))

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = cost_functions["Target One Spin"]

    # Finite difference
    ΔM = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_gradients["Target One Spin"](iso_test)[2:end,:]

    @test round.(fd_M, digits=5) == round.(true_grad, digits=5)
end