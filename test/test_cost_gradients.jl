using Test

@safetestset "Bx: True Control Field Gradient vs Finite Difference" begin
    using GrapeMR
    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = [0.5];
    T2 = [0.25];
    
    # Initial RF field
    N   = 500;
    t_c = 1.0;
    B1  = 10.0;
    B1x, B1y = B1*ones(Float64, 1, N), B1*zeros(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

 ###
    N   = 500;
    αx  = π/2;
    t_c = 0.2;
 
    time = range(0.0, t_c, N);
    t    = time .- t_c/2;
    rotx = rad2deg(αx)/360;
    flip_x = rotx/diff(t)[1];

    BW_Hz = 500.0;
    x     = BW_Hz.*t;
    B0    = zeros(1, N);
    B1x   = ((flip_x.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
    B1y   = zeros(1, N);
    
    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)
    
    mag_test  = forward_propagation(field_test, spin_test[1])
    iso_test  = Magnetization(mag_test, spin_test[1])
    cost      = "Target One Spin"
    cost_func = cost_functions[cost]

    # Finite difference
    Δcf   = 1e-3;
    fd_cf = finite_difference_field(cost_func, field_test, spin_test[1], Δcf, "B1x")

    # True Gradient
    true_grad = get_gradient(field_test, iso_test, Ix, cost)
    
    plot(fd_cf')
    plot!(true_grad')

    @test round.(fd_cf, digits=2) == round.(true_grad, digits=2)

end




@safetestset "By: True Control Field Gradient vs Finite Difference" begin
    using GrapeMR

    # Defining parameters
    M0 = [0.0; 0.0; 1.0];
    T1 = [0.5];
    T2 = [0.25];
    
    # Initial RF field
    N   = 500;
    t_c = 1.0;
    B1  = 10.0;
    B1x, B1y = B1*zeros(Float64, 1, N), B1*ones(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

    ###
    N   = 500;
    αy  = π/2;
    t_c = 0.2;

    time = range(0.0, t_c, N);
    t    = time .- t_c/2;
    roty = rad2deg(αy)/360;
    flip_x = rotx/diff(t)[1];

    BW_Hz = 500.0;
    y     = BW_Hz.*t;
    B0    = zeros(1, N);
    B1x   = zeros(1, N);
    B1y   = ((flip_y.*sinc.(y))./2π)'; # sinc(x) = sin(πx)/(πx)
    
    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

    mag_test  = forward_propagation(field_test, spin_test[1])
    iso_test  = Magnetization(mag_test, spin_test[1])
    cost      = "Target One Spin"
    cost_func = cost_functions[cost]

    # Finite difference
    Δcf   = 1e-3;
    fd_cf = finite_difference_field(cost_func, field_test, spin_test[1], Δcf, "B1y")

    # True Gradient
    true_grad = get_gradient(field_test, iso_test, Iy, cost)
    
    plot(fd_cf')
    plot!(true_grad')

    @test round.(fd_cf, digits=2) == round.(true_grad, digits=2)

end

