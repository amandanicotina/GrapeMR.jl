using Test

@safetestset "Exact ODE solution - Relaxation" begin
    using GrapeMR
    M0 = [1.0; 0.0; 0.0];
    T1 = 0.5;
    T2 = 0.25;
    
    ##### INITIAL RF FIELD #####
    N   = 500;
    α   = π/2;
    t_c = 0.5;
    B1  = α/t_c;
    B1x, B1y = B1*ones(Float64, 1, N), B1*ones(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

    #init_field_test = ControlFields(γ_¹H*B1x, B1y, 0.0, 0.0, t_c, B0, ΔB0)
    
    # Solution Bloch Methods - ODE
    mag = forward_propagation(field_test, spin_test)
    Mx_ODE = mag[2, end]
    My_ODE = mag[3, end]
    Mz_ODE = mag[4, end]

    # Solution for Mz
    time = range(0.0, field_test.t_control, length=length(field_test.B1x)+1)
    Mz_sol = 1.0 .- (1.0 .- spin_test.M_init[3])*exp.(-(2π*spin_test.Γ1).*time)
    #Mz_sol = 1.0 - (1.0 - spin_test.M_init[3])*exp(-field_test.t_control*spin_test.T1)
plot(Mz_sol)
plot!(mag[4,:])

    # Solution for Mxy
    Mxy = √(spin_test.M_init[1]^2 + spin_test.M_init[2]^2)
    #Mxy_sol = Mxy*exp(-t_c/spin_test.T2)
    Mxy_sol = Mxy*exp.(-(2π*spin_test.Γ2).*time)

plot(abs.(Mxy_sol))
plot!(mag[2,:])
plot(angle.(Mxy_sol))
plot!(mag[3,:])

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