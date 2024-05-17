using Test

@safetestset "Test Cost Function: Euclidean Norm" begin
    using GrapeMR
    # Parameters
    N   = 400;
    t_c = 0.1; #[s]

    # Spin
    M0  = [0.0; 0.0; 1.0];
    T1  = [0.6];
    T2  = [0.2];
    B0  = [0.0];
    target = ["max"];
    label  = ["T1-500ms"];

    # RFs
    ΔB1, Bz = 1.0, zeros(1,N);
    B1x = 8*initial_field_spline(N, t_c)'; # rand(1, N); #
    B1y = rand(1,N); #5*initial_field_spline(N, t_c)'; # 

    # Spin and RF objects
    control_field_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_test         = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag_test = forward_propagation(control_field_test, spins_test[1])
    dyn_test = Magnetization(mag_test)
    iso_test = Isochromat(dyn_test, spins_test[1])

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = :euclidean_norm
    cost_val  = GrapeMR.cost_function(iso_test, cost_func)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_function_gradient(iso_test, cost_func)[2:end,:] 

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end


@safetestset "Test Cost Function: Target One Spin" begin
    using GrapeMR
    # Parameters
    N   = 500;
    t_c = 0.1; #[s]

    # Spin
    M0  = [0.0; 0.0; 1.0];
    T1  = [0.5];
    T2  = [0.25];
    B0  = [0.0];
    target = ["max"];
    label  = ["T1-$(Int(T1[1]*1e3))ms"]

    # RFs
    ΔB1, Bz = 1.0, zeros(1,N);
    B1x = initial_field_spline(N, t_c)'; # rand(1, N); #
    B1y = zeros(1,N); #5*initial_field_spline(N, t_c)'; # 

    # Spin and RF objects
    control_field_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_test         = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag_test = forward_propagation(control_field_test, spins_test[1])
    dyn_test = Magnetization(mag_test)
    iso_test = Isochromat(dyn_test, spins_test[1])

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = :target_one_spin
    cost_val  = GrapeMR.cost_function(iso_test, cost_func)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_function_gradient(iso_test, cost_func)[2:end,:] 

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end

@safetestset "Test Cost Function: Saturation Contrast" begin
    using GrapeMR
    # Parameters
    N   = 500;
    t_c = 0.5; #[s]

    # Spin
    M0  = [0.0; 0.0; 1.0];
    T1  = [0.8];
    T2  = [0.6];
    B0  = [0.0];
    target = ["max"];
    label  = ["T1-$(Int(T1[1]*1e3))ms"]

    # RFs
    ΔB1, Bz = 1.0, zeros(1,N);
    B1x = 8*initial_field_spline(N, t_c)'; # rand(1, N); #
    B1y = rand(1,N); #5*initial_field_spline(N, t_c)'; # 

    # Spin and RF objects
    control_field_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_test         = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag_test = forward_propagation(control_field_test, spins_test[1])
    dyn_test = Magnetization(mag_test)
    iso_test = Isochromat(dyn_test, spins_test[1])

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func_sc = :saturation_contrast
    cost_val     = GrapeMR.cost_function(iso_test, cost_func_sc)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func_sc, iso_test, ΔM)

    # True Gradient
    true_grad = cost_function_gradient(iso_test, cost_func_sc)[2:end,:]

    @test round.(fd_M, digits=5) .== round.(true_grad, digits=5)
end