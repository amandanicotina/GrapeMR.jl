using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
T1 = [0.6]
T2 = [0.3]
ΔB1 = [1.0]
B0 = 15.0
label  = ["s1"]
target = ["saturation"]
offset = -B0:1:B0
spins = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
max_iter = 10
poly_start = 0.5
poly_degree = 1
opt_params = OptimizationParams(poly_start, poly_degree, max_iter)

# Grape Parameters 
N = 2000
cost = GrapeMR.euclidean_norm
fields2opt = Dict("B1x" => true, "B1y" => true, "Bz" => false)
grape_params = GrapeParams(N, cost, fields2opt)

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 5.0
control_time = 0.5
control_field = spline_RF(N, control_time, B1ref)

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins)

plot_transverse_magnetization(grape_output.isochromats)
plot_magnetization_2D(grape_output.isochromats)
plot_magnetization_3D(grape_output.isochromats)

plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)
plot_magnetization_time(grape_output.isochromats[2], grape_output.control_field.t_control)
plot_control_fields(grape_output.control_field;unit= "Hz")
plot_control_fields_phase_shift(grape_output.control_field)
plot_transverse_time(grape_output.isochromats, 0.5)

# Cost Function Analysis
GrapeMR.run_cost_analysis(grape_output.control_field, spins[2], 50.0, 20, grape_params.cost_function)

# Save Optimization Data
folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
save_grape_data(grape_output; folder_path)
