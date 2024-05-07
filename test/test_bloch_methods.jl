using Test

@safetestset "Exact ODE solution - Relaxation" begin
    using GrapeMR
    M0 = [0.0; 0.0; 0.074];
    T1 = 0.5;
    T2 = 0.25;
    
    ##### INITIAL RF FIELD #####
    N   = 250;
    t_c = 0.5;
    B1  = 1.0;
    B1x, B1y = B1*zeros(Float64, 1, N), B1*zeros(Float64, 1, N)
    B0, ΔB0  = [100.0, 150.0], zeros(1, N);

    (spins, field_test) = normalization(M0, T1, T2, B0, target, label, t_c, B1x, B1y, B1, Bz)
    
    # Solution Bloch Methods - ODE
    spin_test = spins[2]
    mag_test  = forward_propagation(field_test, spin_test)
    Mx_ODE  = mag_test[2, :]
    My_ODE  = mag_test[3, :]
    Mz_ODE  = mag_test[4, :]

    # Solution for Mz
    time = range(0.0, field_test.t_control, length=length(field_test.B1x)+1)
    Mz_sol = 1.0 .- (1.0 .- spin_test.M_init[3])*exp.(-(2π*spin_test.T1).*time)

    # Solution for Mxy
    Mxy_ini = Mxy_ODE[1]
    Mx_sol = Mx_ODE*exp.(-(2π*spin_test.T2).*time).*exp.(-(2π*im*spin_test.B0inho).*time)
    My_sol = My_ODE*exp.(-(2π*spin_test.T2).*time).*exp.(-(2π*im*spin_test.B0inho).*time)

    # Plots
    plot(time, Mz_ODE)
    plot!(time, Mz_sol)

    plot(Mxy_ODE)
    plot!(Mxy_sol)

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

@safetestset "Bloch Equation units vs normalized" begin
    # Units
    t_c = 1.0;
    N   = 500;
    B0  = [0.0];
    Bz  = zeros(1,N);
    B1x = [ones(1, 250) zeros(1, 250)];#initial_field_spline(N, t_c)'; 
    B1y = zeros(1,N);#initial_field_spline(N, t_c)';
    spin = Spin([0.0, 0.0, 1.0], 0.5, 0.25, 0.0, 0.0, 0.0, "test", "test", 1.0)
    control_field = ControlField(B1x, B1y, 1.0, Bz, t_c)
    function bloch_mat(B1x::Float64, B1y::Float64, Bz::Float64, T1::Float64, T2::Float64)
        bloch_matrix = 
            [0.0   0.0   0.0   0.0;
             0.0  -1/T2  Bz   -B1y;
             0.0  -Bz   -1/T2  B1x;
             1/T1  B1y  -B1x  -1/T1] 
        
        return bloch_matrix
    end   
    function forward_prop(cf::ControlField, s::Spin)
        Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
        M       = zeros(Float64, 4, length(cf.B1x)+1)
        M[:, 1] = [1.0; s.M_init[1]; s.M_init[2]; s.M_init[3]];
        
        B0 = s.B0inho
        B1 = s.B1inho
        Bz = cf.Bz .+ B0
        Bx = B1*cf.B1x
        By = B1*cf.B1y
    
        for (i, Δt) ∈ enumerate(diff(Δt_arr))
            b_m = 2π*bloch_mat(Bx[i], By[i], Bz[i], s.T1, s.T2)
            M[:, i+1] = expv(Δt, b_m, M[:, i])
        end

        return M    
    end
    mag_unit  = forward_prop(control_field, spin)
    mag_norm = forward_propagation(control_field, spin)
    time = range(0.0, t_c, length=N+1)
    using Plots
    plot(time, mag_unit[2:end, :]')
    plot(time, mag_norm[2:end, :]')
end