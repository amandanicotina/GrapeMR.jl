using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = [0.0] 
T1 = [0.6, 0.1]
T2 = [0.3, 0.05]
label = ["s1", "s2"]
target = ["min", "max"]
spins = Spin(M0, T1, T2, offsets, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.saturation_contrast, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
random_opt = @time random_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), range(500, 1000, step = 250)) 
bohb_opt   = @time bohb_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 1000) 

## HPO Plots ##
plot_hyperopt_history(random_opt; title = "Random Sampler") 
plot_cost_grape_runs(random_opt; plotlog=false, title = "Random Sampler") 
plot_cost_hyperparam(random_opt; grid_plot=true, title = "Random Sampler")

plot_hyperopt_history(bohb_opt; title = "BOHB Sampler") 
plot_cost_grape_runs(bohb_opt; plotlog=true, title = "BOHB Sampler") 
plot_cost_hyperparam(bohb_opt; grid_plot=true, title = "BOHB Sampler") 
## HPO Plots ##

Tc, poly_start, poly_degree, max_iter = bohb_opt.minimizer;
opt_params = OptimizationParams(poly_start, poly_degree, max_iter);

# Parameters 
params = Parameters(grape_params, opt_params);

# Initial RF Pulse
B1ref = 1.0;
control_field = spline_RF(grape_params.N, Tc, B1ref);

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins);

# Plots
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)

# Cost Function Analysis
GrapeMR.run_cost_analysis(grape_output.control_field, spins[1], 80.0, 20, grape_output.params.grape_params.cost_function)

# # Save Optimization Data
# save_grape_data(grape_output; folder_path)
# save_hyperopt_data(random_opt; folder_path)