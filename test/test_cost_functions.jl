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
    dyn_test = GrapeMR.Magnetization(mag_test)
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
    dyn_test = GrapeMR.Magnetization(mag_test)
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
    control_sc_test = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins_sc_test   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag_sc_test = forward_propagation(control_sc_test, spins_sc_test[1])
    dyn_sc_test = GrapeMR.Magnetization(mag_sc_test)
    iso_sc_test = Isochromat(dyn_sc_test, spins_test[1])

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

# New cost function test: steady State

using GrapeMR
using BlochSim
using Plots

# Parameters bSSFP
N  = 401                 # Number of points
α  = π / 8               # Flip angle in radians
Δϕ = π                  # Phase cycling
TR = 5e-3               # Repetition time in seconds

# Parameters RF
t_c = 0.5;
ΔB1, Bz  = 1.0, zeros(1, N)
B1x, B1y = zeros(1, N), zeros(1, N)

# Spin parameters
M0 = [0.0, 0.0, 1.0]    # Initial magnetization vector
T1 = [1.0]              # Longitudinal relaxation time
T2 = [0.8]              # Transverse relaxation time
B0 = range(-2/TR, stop=2/TR, length=N) |> collect # B0 field range

# Target and label for simulation
target = ["max"]
label  = ["T1-$(Int(T1[1]*1e3))ms"]

# Spins and ControlField objects
spins_ss_test   = GrapeMR.SteadyState(M0, T1, T2, B0, ΔB1, target, label, α, Δϕ, TR, TR/2)
control_ss_test = ControlField(B1x, B1y, 1.0, Bz, t_c)

mag_ss_test = forward_propagation(control_ss_test, spins_ss_test[1])
dyn_ss_test = GrapeMR.Magnetization(mag_ss_test)
iso_ss_test = Isochromat(dyn_ss_test, spins_ss_test[1]) 

cost_func_ss = :target_steady_state
cost_val_ss  = GrapeMR.cost_function(iso_ss_test, cost_func_ss)

# Finite difference
ΔM = 1e-10
fd_M_ss = finite_difference_cost(cost_func_ss, iso_ss_test, ΔM)

# True gradient
true_grad_ss = cost_function_gradient(iso_ss_test, cost_func_ss)[2:end, :]

# using ForwardDiff
# ss = steady_state_matrix(spins_ss_test[1])
# Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
# ss_vec  = [Mx_ss, My_ss, Mz_ss]
# Mx, My, Mz = mag_ss_test[2,end], mag_ss_test[3,end], mag_ss_test[4,end]
# mag_vec = [Mx, My, Mz]
# gradient = ForwardDiff.gradient()