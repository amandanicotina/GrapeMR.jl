using Test

@safetestset "Exact ODE solution - Relaxation" begin
    using GrapeMR
    spin_test = Spins([1.0; 0.0; 0.0], 0.6, 0.3, 0.0, "nothing")
    #N   = 100
    #t_c = 0.5
    #B1x, B1y = zeros(Float64, 1, N), zeros(Float64, 1, N)
    #B0, ΔB0 = zeros(1, N), zeros(1, N);

    N    = 2000;
    α    = π/2;
    t_c  = 0.1;
    B1x_max, B1y_max = 2π/t_c, 2π/t_c;

    time = range(0.0, t_c, N) 
    t = time .- t_c/2;
    flip_rads = α/diff(t)[1]
    BW_Hz = 500;
    x = 2*BW_Hz.*t/2
    y = BW_Hz.*t/2
    B1x_arr = flip_rads*sinc.(x) # sinc(x) = sin(πx)/(πx)
    B1y_arr = flip_rads*sinc.(y) #  

    B0  = zeros(1, N);
    ΔB0 = zeros(1, N);
    field_test = ControlFields((B1x_arr./maximum(B1x_arr))', (B1y_arr./maximum(B1x_arr))', 
                    B1x_max, B1y_max, t_c, B0, ΔB0); 

    #init_field_test = ControlFields(γ_¹H*B1x, B1y, 0.0, 0.0, t_c, B0, ΔB0)
    
    # Solution Bloch Methods - ODE
    mag = forward_propagation(field_test, spin_test)
    Mx_ODE = mag[2, end]
    My_ODE = mag[3, end]
    Mz_ODE = mag[4, end]

    # Solution for Mz
    # Mz_sol = 1.0 .- (1.0 .- spin_test.M_init[3])*exp.(-(2π/spin_test.T1).*time)
    Mz_sol = 1.0 - (1.0 - spin_test.M_init[3])*exp(-2π*t_c/spin_test.T1)

    # Solution for Mxy
    Mxy = √(spin_test.M_init[1]^2 + spin_test.M_init[2]^2)
    Mxy_sol = Mxy*exp(-2π*t_c/spin_test.T2)
    #Mxy_sol = Mxy*exp.(-(2π/spin_test.T2).*time)

    Mx_sol = abs(Mxy_sol)
    My_sol = 0.0

    @test round(Mx_ODE, digits=5) == round(Mx_sol, digits=5)
    @test round(My_ODE, digits=5) == round(My_sol, digits=5)
    @test round(Mz_ODE, digits=5) == round(Mz_sol, digits=5)
end

@safetestset "Exact ODE solution - Rotation" begin
    using GrapeMR
    spin_test = Spins([0.0; 0.0; 1.0], 1e8, 1e8, 0.0, "nothing")
    M_ini_sol = [spin_test.M_init[1]; spin_test.M_init[2]; spin_test.M_init[3]]

    # Flip angle and RF field
    N   = 100 
    α   = π/3;
    t_c = 1e-6;
    B1  = α/(γ_¹H*t_c);
    B1x, B1y = B1*ones(Float64, 1, N), zeros(Float64, 1, N)
    B0, ΔB0 = zeros(1, N), zeros(1, N);
    init_field_test = ControlFields(γ_¹H*B1x, B1y, 0.0, 0.0, t_c, B0, ΔB0)

    # Solution Bloch Methods - ODE
    mag = forward_propagation(init_field_test, spin_test)
    Mx_ODE = mag[2, end]
    My_ODE = mag[3, end]
    Mz_ODE = mag[4, end]

    # Solution for Mz
    Rx = [1.0  0.0     0.0;
          0.0  cos(α)  sin(α);
          0.0 -sin(α)  cos(α)]
    M_sol = Rx*M_ini_sol

    @test round(Mx_ODE, digits=4) == round(M_sol[1], digits=4)
    @test round(My_ODE, digits=4) == round(M_sol[2], digits=4)
    @test round(Mz_ODE, digits=4) == round(M_sol[3], digits=4)
end