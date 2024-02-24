using Test

@safetestset "Test Cost Function: Euclidean Norm" begin
    using GrapeMR
    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = 0.5;
    T2 = 0.25;
    
    ##### INITIAL RF FIELD #####
    N   = 500;
    t_c = 1.0;
    B1  = 10.0;
    B1x, B1y = B1*ones(Float64, 1, N), B1*ones(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)
    
    mag_test = forward_propagation(field_test, spin_test)
    iso_test = Magnetization((mag_test,), (spin_test,))

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = cost_functions["Euclidean Norm"]
    cost_func(iso_test)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_gradients["Euclidean Norm"](iso_test)[2:end,:]

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end


@safetestset "Test Cost Function: Target One Spin" begin
    using GrapeMR
    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = 0.5;
    T2 = 0.25;
    
    ##### INITIAL RF FIELD #####
    N   = 500;
    t_c = 1.0;
    B1  = 10.0;
    B1x, B1y = B1*ones(Float64, 1, N), B1*ones(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)
    
    mag_test = forward_propagation(field_test, spin_test)
    iso_test = Magnetization((mag_test,), (spin_test,))

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = cost_functions["Target One Spin"]
    cost_func(iso_test)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_gradients["Target One Spin"](iso_test)[2:end,:]

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end