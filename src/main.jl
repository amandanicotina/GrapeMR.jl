using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = -15:1:15
# offsets = [0.0]
T1 = [0.6, 0.1]
T2 = [0.3, 0.05]
label = ["s1", "s2"]
target = ["max", "min"]
spins = Spin(M0, T1, T2, offsets, ΔB1, target, label)

# Optimization Parameters
poly_start = 0.5
poly_degree = 1
max_iter = 25000
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.saturation_contrast, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
Tc = 0.1
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins)

# Plots
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)

# Cost Function Analysis
GrapeMR.run_cost_analysis(grape_output.control_field, spins[end], 80.0, 20, grape_params.cost_function)

# Save Optimization Data
folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
save_grape_data(grape_output; folder_path)