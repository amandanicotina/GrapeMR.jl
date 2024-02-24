using Test

@safetestset "Exact ODE solution - Relaxation" begin
    using GrapeMR
    M0 = [0.5; 0.5; 0.0];
    T1 = 0.5;
    T2 = 0.25;
    
    ##### INITIAL RF FIELD #####
    N   = 500;
    t_c = 1.0;
    B1x, B1y = B1*zeros(Float64, 1, N), B1*zeros(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

    #init_field_test = ControlFields(γ_¹H*B1x, B1y, 0.0, 0.0, t_c, B0, ΔB0)
    
    # Solution Bloch Methods - ODE
    mag = forward_propagation(field_test, spin_test)
    Mx_ODE  = mag[2, :]
    My_ODE  = mag[3, :]
    Mz_ODE  = mag[4, :]
    Mxy_ODE = sqrt.((Mx_ODE).^2 .+ (My_ODE).^2)

    # Solution for Mz
    time = range(0.0, field_test.t_control, length=length(field_test.B1x)+1)
    Mz_sol = 1.0 .- (1.0 .- spin_test.M_init[3])*exp.(-(2π*spin_test.T1).*time)

    # Solution for Mxy
    Mxy = √(spin_test.M_init[1]^2 + spin_test.M_init[2]^2)
    Mxy_sol = Mxy*exp.(-(2π*spin_test.T2).*time)

    @test round.(Mxy_ODE, digits=5).== round.(Mxy_sol, digits=5)
    @test round.(Mz_ODE, digits=5) .== round.(Mz_sol, digits=5)
    #@test round(My_ODE, digits=5) == round(My_sol, digits=5)
end

@safetestset "Exact ODE solution - Rotation" begin
    using GrapeMR
    M0 = [0.0; 0.0; 1.0];
    T1 = 1e8;
    T2 = 1e8;
    
    # Flip angle and RF field
    N   = 500 
    α   = π/2;
    t_c = 1e-6;
    rot = rad2deg(α)/360
    B1  = rot/t_c;
    B1x, B1y = B1*ones(Float64, 1, N), B1*zeros(Float64, 1, N)
    B0, ΔB0  = zeros(1, N), zeros(1, N);

    (spin_test, init_field_test) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

    # Solution Bloch Methods - ODE
    mag = forward_propagation(init_field_test, spin_test)
    Mx_ODE = mag[2, end]
    My_ODE = mag[3, end]
    Mz_ODE = mag[4, end]

    # Solution for Mz
    Rx = [1.0  0.0     0.0;
          0.0  cos(α)  -sin(α);
          0.0  sin(α)  cos(α)]

    Ry = [cos(α)  0.0   sin(α);
          0.0     1.0   0.0;
          -sin(α) 0.0   cos(α)]

    M_sol = Rx*M0

    @test round(Mx_ODE, digits=4) == round(M_sol[1], digits=4)
    @test round(My_ODE, digits=4) == round(M_sol[2], digits=4)
    @test round(Mz_ODE, digits=4) == round(M_sol[3], digits=4)
end