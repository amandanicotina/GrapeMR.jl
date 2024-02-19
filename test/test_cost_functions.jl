using Test

@safetestset "True Cost Function Gradient vs Finite Difference" begin
    using GrapeMR
    # Defining parameters
    #N   = 1000;
    #α   = π/2;
    #t_c = 1.0;
    #B1  = α/(γ_¹H*t_c);
    #B1x, B1y = B1*ones(1, N), B1*zeros(1, N);   
    #B0, ΔB0  = zeros(1, N), zeros(1, N);


    #field_test = ControlFields(B1x, B1y, 
    #B1x_max, B1y_max, t_c, B0, ΔB0); 
    N    = 500;
    α    = π/2;
    t_c  = 0.25;
    B1x_max, B1y_max = 2π/t_c, 2π/t_c;
                    
    time = range(0.0, t_c, N) 
    t = time .- t_c/2;
    flip_rads = α/diff(t)[1]
    BW_Hz = 500;
    x = 2*BW_Hz.*t/2
    y = BW_Hz.*t/2
    B1x_arr = flip_rads*sinc.(x); # sinc(x) = sin(πx)/(πx)
    B1y_arr = flip_rads*sinc.(y); #;  
                    
    B0  = zeros(1, N);
    ΔB0 = zeros(1, N);
                    
    field_test = ControlFields((B1x_arr./maximum(B1x_arr))', (B1y_arr./maximum(B1y_arr))', B1x_max, B1y_max, t_c, B0, ΔB0); 
                    
    spin_test  = Spins([0.0; 0.0; 1.0], 0.5, 0.3, 0.0, "max")
    mag_test   = forward_propagation(field_test, spin_test)
    iso_test   = Magnetization((mag_test,), (spin_test,))

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];
    M_tar = [0.5, 0.5, 0.0]

    cost_func = cost_functions["Euclidean Norm"]
    cost_func(iso_test)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    # True Gradient
    true_grad = cost_gradients["Euclidean Norm"](iso_test)[2:end,:]

    @test round.(fd_M, digits=5) == round.(true_grad, digits=5)
end