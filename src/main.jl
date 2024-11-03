using GrapeMR

# Spin Parameters
<<<<<<< Updated upstream
M0      = [0.0, 0.0, 1.0]
ΔB1     = [1.0]
offsets = -15:1:15
T1 = [0.6, 0.1]
T2 = [0.3, 0.05]
label  = ["s1", "s2"]
target = ["max", "min"]
spins  = Spin(M0, T1, T2, offsets, ΔB1, target, label)

# Optimization Parameters
Tc = 1.0
poly_start  = 0.5
poly_degree = 1
max_iter    = 6000
opt_params  = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.saturation_contrast, [true true false])
=======
M0     = [0.0, 0.0, 1.0]
ΔB1    = [1.0]
# offset = -1:0.1:1
offset = 1:1:1
T1 = [1.35]
T2 = [0.05]
label  = ["Blood"]
target = ["-"]
spins  = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
Tc = 0.1
poly_start  = 0.1
poly_degree = 1.0
max_iter    = 1
opt_params  = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.euclidean_norm, Dict("B1x" => true, "B1y" => true, "Bz" => false))
>>>>>>> Stashed changes

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 5.0
control_field = hard_RF(grape_params.N, Tc, B1ref) 

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins)

plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)

# spin = grape_output.isochromats[1].spin
# @time run_cost_analysis(grape_output.control_field, spin, 150.0, 50, grape_params.cost_function) 

# # folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
# # save_grape_data(grape_output; folder_path)

# using GrapeMR
# load_folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/2024-11-02/s1_saturation_contrast_min_15Hz"
# grape_output1 = load_grape_data(load_folder_path);
# plot_magnetization_control_field(grape_output1.control_field, grape_output1.isochromats)
# plot_cost_values(grape_output1.cost_values, grape_output1.params.grape_params)

# spin = grape_output1.isochromats[1].spin;
# @time run_cost_analysis(grape_output1.control_field, spin, 500.0, 50, grape_output1.params.grape_params.cost_function)


# # Spin System
# M0 = [0.0, 0.0, 1.0]
# ΔB1 = [1.0]
# B0 = 0.0
# target = ["max", "min"]
# label  = ["s1", "s2"]

# T1 = [0.5, 0.2]
# T2 = [0.1, 0.05]

# α, Δϕ, TR, TE = π/3, 2π, 5e-3, 2.5e-3
# spins = GrapeMR.SteadyState(M0, T1, T2, B0, ΔB1, target, label, α, Δϕ, TR, TE)

# # Optimization Parameters
# Tc = 0.5
# poly_start = 0.1
# poly_degree = 1.0
# max_iter = 1000
# opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# # Grape Parameters 
# grape_params = GrapeParams(1500, GrapeMR.saturation_contrast_steady_state, [true true false])

# # Parameters 
# params = Parameters(grape_params, opt_params)

# # Initial RF Pulse
# B1ref = 1.0
# control_field = spline_RF(grape_params.N, Tc, B1ref) 

# # Run Optimization
# grape_output = @time grape(params, control_field, spins); 

# # Plots
# plot_cost_values(grape_output.cost_values, grape_params)
# plot_magnetization_target(grape_output.isochromats)
# plot_Mtrans_offset_ss(grape_output.isochromats)
# plot_magnetization_2D(grape_output.isochromats)


