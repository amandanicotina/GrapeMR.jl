using GrapeMR, BenchmarkTools

# Spin Parameters
M0     = [0.0, 0.0, 1.0]
ΔB1    = [1.0]
offset = [0.0]
T1 = [1.35]
T2 = [0.05]
label  = ["Blood"]
target = ["-"]
spins  = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
Tc = 0.1
poly_start  = 0.1
poly_degree = 1.0
max_iter    = 1500.0
opt_params  = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, :euclidean_norm, [true true false])

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref) 

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins); 
# @btime grape(params, control_field, spins)
# @benchmark grape(params, control_field, spins)

# Plots
plot_cost_values(grape_output.cost_values, grape_params)
plot_magnetization_2D(grape_output.isochromats)
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)




### Steady State ###
# plot_ss_offset_profile(spins)
# flip_angles = [π/4, π/3, π/2, 3π/4, π]
# p = plot_ss_flip_angle_profile(spins, flip_angles)
# Spin System
# M0 = [0.0, 0.0, 1.0]
# ΔB1 = [1.0]
# B0 = 50.0
# target = ["max"]
# label  = ["s1"]
# T1 = [0.5]
# T2 = [0.1]
# offset = collect(-B0:10:B0)
# α, Δϕ, TR, TE = π/3, 2π, 5e-3, 2.5e-3
# spins = GrapeMR.SteadyState(M0, T1, T2, offset, ΔB1, target, label, α, Δϕ, TR, TE)

# plot_ss_offset_profile(spins)

# plot_magnetization_target(grape_output.isochromats)

# plot_magnetization_targetB0(grape_output.isochromats)

# plot_Mtrans_offset_ss(grape_output.isochromats)

# plot_magnetization_3D(grape_output.isochromats)

# plot_magnetization_target_3D(grape_output.isochromats[1])
