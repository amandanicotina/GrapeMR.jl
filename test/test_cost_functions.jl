using Test

@safetestset "Test Cost Function: Euclidean Norm" begin
    using GrapeMR
    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = [0.5];
    T2 = [0.2];
    target = ["min"];
    
    ##### INITIAL RF FIELD #####
    N   = 500;
    t_c = 1.0;
    B1  = 10.0;
    B1x, B1y = B1*ones(Float64, 1, N), B1*ones(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, T1, T2, target, t_c, B1x, B1y, B0);
    
    mag_test = forward_propagation(field_test, spin_test[1])
    dyn_test = Magnetization(mag_test)
    iso_test = Isochromat(dyn_test, spin_test[1])

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = euclidean_norm(iso_test)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(euclidean_norm, iso_test, ΔM)

    # True Gradient
    true_grad = grad_euclidean_norm(iso_test)[2:end,:] 

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end


@safetestset "Test Cost Function: Target One Spin" begin
    using GrapeMR
    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = [0.5];
    T2 = [0.2];
    target = ["min"];
    
    ##### INITIAL RF FIELD #####
    N   = 500;
    t_c = 1.0;
    B1  = 10.0;
    B1x, B1y = B1*ones(Float64, 1, N), B1*ones(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, T1, T2, target, t_c, B1x, B1y, B0);
    
    mag_test = forward_propagation(field_test, spin_test[1])
    dyn_test = Magnetization(mag_test)
    iso_test = Isochromat(dyn_test, spin_test[1])

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = target_one_spin(iso_test)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(target_one_spin, iso_test, ΔM)

    # True Gradient
    true_grad = grad_target_one_spin(iso_test)[2:end,:]

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end

@safetestset "Test Cost Function: Saturation Contrast" begin
    using GrapeMR
    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = [0.1];
    T2 = [0.01];
    target = ["max"];
    
    ##### INITIAL RF FIELD #####
    t_c = 0.5;
    N   = 600;
    B0  = [0.0];
    #B0  = [-150.0, -100.0, -50, 0.0, 50.0, 100.0, 150.0]; # [Hz]
    Bz  = zeros(1,N);
    B1_max = 500; # [Hz]
    B1x = B1_max*initial_field_spline(N, t_c)'; 
    B1y = B1_max*initial_field_spline(N, t_c)';

    (spins, field_init) = normalization(M0, T1, T2, B0, target, label, t_c, B1x, B1y, B1, Bz);

    spin_test = spins[1]
    mag_test  = forward_propagation(field_test, spin_test)
    dyn_test  = Magnetization(mag_test)
    iso_test  = Isochromat(dyn_test, spin_test)

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = cost_saturation_contrast(iso_test)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_saturation_contrast, iso_test, ΔM)

    # True Gradient
    true_grad = grad_saturation_contrast(iso_test)[2:end,:]

    @test round.(fd_M, digits=5) == round.(true_grad, digits=5)
end