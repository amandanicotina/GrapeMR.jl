using Test

@safetestset "Exact ODE solution - Relaxation" begin
    using GrapeMR
    spin_test = Spins([1.0; 0.0; 0.0], 0.6, 0.3, 0.0, "nothing")
    N   = 100
    t_c = 0.5
    B1x, B1y = zeros(Float64, 1, N), zeros(Float64, 1, N)
    init_field_test = InitialControlFields(N, B1x, 0.0, B1y, 0.0, t_c, 0.0, 0.0)
    
    # Solution Bloch Methods - ODE
    mag = magnetization_ODE(init_field_test, spin_test)
    Mx_ODE = mag[2, end]
    My_ODE = mag[3, end]
    Mz_ODE = mag[4, end]

    # Solution for Mz
    Mz_sol = 1.0 - (1.0 - spin_test.M_init[3])*exp(-t_c/spin_test.T1)

    # Solution for Mxy
    Mxy = √(spin_test.M_init[1]^2 + spin_test.M_init[2]^2)
    Mxy_sol = Mxy*exp(-t_c/spin_test.T2)

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
    t_c = 0.001;
    B1  = α/(γ_¹H*t_c);
    B1x, B1y = B1*ones(Float64, 1, N), zeros(Float64, 1, N)
    init_field_test = InitialControlFields(N, B1x, 0.0, B1y, 0.0, t_c, 0.0, 0.0)

    # Solution Bloch Methods - ODE
    mag = magnetization_ODE(init_field_test, spin_test)
    Mx_ODE = mag[2, end]
    My_ODE = mag[3, end]
    Mz_ODE = mag[4, end]

    # Solution for Mz
    Rx = [1.0  0.0      0.0;
          0.0  cos(α)   -sin(α);
          0.0  sin(α)   cos(α)]
    M_sol = Rx*M_ini_sol

    @test round(Mx_ODE, digits=5) == round(M_sol[1], digits=5)
    @test round(My_ODE, digits=5) == round(M_sol[2], digits=5)
    @test round(Mz_ODE, digits=5) == round(M_sol[3], digits=5)
end