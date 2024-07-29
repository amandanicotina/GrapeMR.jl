using Test

@safetestset "Test Cost Function: Euclidean Norm" begin
    using GrapeMR
    # Parameters
    N   = 400;
    t_c = 0.5; #[s]

    # Spin
    M0  = [0.0; 0.0; 1.0];
    T1  = [0.6];
    T2  = [0.2];
    B0  = [0.0];
    target = ["max"];
    label  = ["T1-500ms"];

    # RFs
    ΔB1, Bz = 1.0, zeros(1,N);
    B1x = 8*spline_RF(N, t_c)'; # rand(1, N); #
    B1y = rand(1,N); #5*initial_field_spline(N, t_c)'; # 

    # Spin and RF objects
    control_field_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_test         = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag_test = forward_propagation(control_field_test, spins_test[1])
    dyn_test = GrapeMR.Magnetization(mag_test)
    iso_test = Isochromat(dyn_test, spins_test[1])

    Mx = mag_test[2,:];
    My = mag_test[3,:];
    Mz = mag_test[4,:];

    cost_func = :euclidean_norm
    (cost_val, cost_grad)  = GrapeMR.cost_function(iso_test, cost_func)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(cost_func, iso_test, ΔM)

    @test round.(fd_M, digits=5) .== round.(cost_grad[2:end,:] , digits=5)
end


@safetestset "Test Cost Function: Target One Spin" begin
    using GrapeMR
    # Parameters
    N   = 1000;
    t_c = 0.5; #[s]
    grape_params = GrapeParams(1000, :target_one_spin, [true true false])

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8];
    T2  = [0.6];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = ["min"];
    label  = ["T1-$(Int(T1[1]*1e3))ms"]

    # RFs
    B1ref = 1.0
    B1x = spline_RF(grape_params.N, t_c)'
    B1y = spline_RF(grape_params.N, t_c)'
    Bz  = zeros(1, grape_params.N)
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)

    # Spin and RF objects
    control_sc_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_sc_test   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)
    
    mag_sc_test = forward_propagation(control_sc_test, spins_sc_test[1])
    dyn_sc_test = GrapeMR.Magnetization(mag_sc_test)
    iso_sc_test = Isochromat(dyn_sc_test, spins_sc_test[1])

    Mx = mag_sc_test[2,:];
    My = mag_sc_test[3,:];
    Mz = mag_sc_test[4,:];

    (cost_val, cost_grad) = GrapeMR.cost_function(iso_sc_test, grape_params.cost_function)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(grape_params.cost_function, iso_sc_test, ΔM)

    @test round.(fd_M, digits=5) .== round.(cost_grad[2:end,:], digits=5)

end

@safetestset "Test Cost Function: Saturation Contrast" begin
    # Parameters
    using GrapeMR
    N   = 1000;
    t_c = 0.5; #[s]
    grape_params = GrapeParams(1000, :saturation_contrast, [true true false])

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8, 0.4];
    T2  = [0.6, 0.3];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = ["min", "max"];
    label  = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

    # RFs
    B1ref = 1.0
    B1x = spline_RF(grape_params.N, t_c)'
    B1y = spline_RF(grape_params.N, t_c)'
    Bz  = zeros(1, grape_params.N)
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)

    # Spin and RF objects
    control_sc_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_sc_test   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)
    
    mag_sc_test = forward_propagation(control_sc_test, spins_sc_test[1])
    dyn_sc_test = GrapeMR.Magnetization(mag_sc_test)
    iso_sc_test = Isochromat(dyn_sc_test, spins_sc_test[1])

    Mx = mag_sc_test[2,:];
    My = mag_sc_test[3,:];
    Mz = mag_sc_test[4,:];

    (cost_val, cost_grad) = GrapeMR.cost_function(iso_sc_test, grape_params.cost_function)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(grape_params.cost_function, iso_sc_test, ΔM)

    @test round.(fd_M, digits=5) .== round.(cost_grad[2:end,:], digits=5)

end 

# @safetestset "Test Cost Function: Saturation Contrast Mx" begin
using GrapeMR
# Parameters
N   = 1000;
t_c = 0.5; #[s]
grape_params = GrapeParams(1000, :target_different_offsets_steady_state, [true true false])

# Spin
M0 = [0.0, 0.0, 1.0] 
T1 = [0.8, 0.08]
T2 = [0.2, 0.04]
label  = ["T1=$(T1[1]*1e3)ms", "T1=$(T1[2]*1e3)ms"]  
target = ["min", "max"]
B0 = 30.0
offset = collect(-B0/2:3:B0/2) 
ΔB1 = [1.0]
α, Δϕ, TR, TE = 2π/9, π, 5e-3, 5e-3/2

# RFs
B1ref = 1.0
B1x = spline_RF(grape_params.N, t_c)'
B1y = spline_RF(grape_params.N, t_c)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)

# Spin and RF objects
control_sc_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
# spins_sc_test   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)
spins_sc_test = GrapeMR.SteadyState(M0, T1, T2, offset, ΔB1, target, label, α, Δϕ, TR, TE)

mag_sc_test = forward_propagation(control_sc_test, spins_sc_test[1])
dyn_sc_test = GrapeMR.Magnetization(mag_sc_test)
iso_sc_test = Isochromat(dyn_sc_test, spins_sc_test[1])

Mx = mag_sc_test[2,:];
My = mag_sc_test[3,:];
Mz = mag_sc_test[4,:];

(cost_val, cost_grad) = GrapeMR.cost_function(iso_sc_test, grape_params.cost_function)

# Finite difference
ΔM   = 1e-10;
fd_M = finite_difference_cost(grape_params.cost_function, iso_sc_test, ΔM)

@test round.(fd_M, digits=5) .== round.(cost_grad[2:end,:], digits=5)
# end