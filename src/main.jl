using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = -15:5:15
T1 = [0.6, 0.1] #[1/31.3436]
T2 = [0.3, 0.05] #[1/37.6471]
label = ["s1", "s2"]
target = ["max", "min"]
spins = Spin(M0, T1, T2, offsets, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.saturation_contrast, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
# random_opt = @time random_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), range(1500, 3000, step = 500)) 
# bohb_opt = @time bohb_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 3000)
# hband_opt = @time hband_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 10)
Tc, poly_start, poly_degree, max_iter = 0.1, 0.7, 1, 30000;
opt_params = OptimizationParams(poly_start, poly_degree, max_iter);

# Parameters 
params = Parameters(grape_params, opt_params);

# Initial RF Pulse
B1ref = 1.0;
control_field = spline_RF(grape_params.N, Tc, B1ref);

# Run Optimization
grape_output6 = @time GrapeMR.grape(params, control_field, spins);
plot_magnetization_control_field(grape_output6.control_field, grape_output6.isochromats)
plot_cost_values(grape_output6.cost_values, grape_output6.params.grape_params)

# Plots
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)
plot_magnetization_3D(grape_output.isochromats)
plot_magnetization_3D(grape_output2.isochromats)
plot_magnetization_control_field(grape_output2.control_field, grape_output2.isochromats)
plot_cost_values(grape_output2.cost_values, grape_output2.params.grape_params)
plot_magnetization_control_field(grape_output3.control_field, grape_output3.isochromats)
plot_cost_values(grape_output3.cost_values, grape_output3.params.grape_params)
plot_magnetization_control_field(grape_output4.control_field, grape_output4.isochromats)
plot_cost_values(grape_output4.cost_values, grape_output4.params.grape_params)


# plot_hyperopt_history(random_opt; title = "Random Sampler") 
# plot_hyperopt_history(bohb_opt; title = "BOHB Sampler") 
# plot_cost_grape_runs(bohb_opt; plotlog = true, title = "BOHB Sampler")
# plot_cost_hyperparam(bohb_opt; grid_plot = true, title = "BOHB Sampler")

# Cost Function Analysis
# GrapeMR.run_cost_analysis(grape_output.control_field, spins[1], 80.0, 20, grape_output.params.grape_params.cost_function)

# Save Optimization Data
folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
save_grape_data(grape_output; folder_path)
save_grape_data(grape_output2; folder_path)
save_grape_data(grape_output3; folder_path)
save_grape_data(grape_output4; folder_path)
save_grape_data(grape_output5; folder_path)
# save_hyperopt_data(random_opt; folder_path)
