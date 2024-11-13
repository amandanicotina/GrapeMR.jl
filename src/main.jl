using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = -50:5:50
# offsets = [0.0]
T1 = [1/31.3436]
T2 = [1/37.6471]
label = ["s1"]
target = ["min"]
no_relax_spins = Spin(M0, [1e8], [1e8], offsets, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.spin_target, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
random_opt_exp = @time random_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), range(1500, 3500, step = 500)) 
random_opt = @time random_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), range(1500, 3000, step = 500)) 
random_opt_inho = @time random_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), range(1000, 5000, step = 500)) 
# 25.836789 seconds (3.91 M allocations: 14.725 GiB, 1.21% gc time)
# LinRange(0.1, 0.5, 10), range(1500, 3000, step = 500) => 2564.183795 seconds (236.47 M allocations: 1.259 TiB, 1.78% gc time, 0.01% compilation time)

bohb_opt_exp = @time bohb_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 3000) 
bohb_opt_inho = @time bohb_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 5000) 
# LinRange(0.1, 0.5, 10), 3000 => 43.065560 seconds (66.77 M allocations: 24.954 GiB, 1.46% gc time)
# LinRange(0.1, 0.5, 10), 3000) => 5635.360429 seconds (712.27 M allocations: 2.741 TiB, 1.35% gc time, 0.07% compilation time)
#8185.092417 seconds (947.28 M allocations: 4.087 TiB, 1.20% gc time, 0.00% compilation time)
gr()
plot_hyperopt_history(random_opt_exp; title = "Random Sampler") 
plot_cost_grape_runs(random_opt_exp; plotlog=false, title = "Random Sampler") 
plot_cost_hyperparam(random_opt_exp; grid_plot=true, title = "Random Sampler")

plot_hyperopt_history(bohb_opt_exp; title = "BOHB Sampler") 
plot_cost_grape_runs(bohb_opt_exp; plotlog=true, title = "BOHB Sampler") 
plot_cost_hyperparam(bohb_opt_exp; grid_plot=true, title = "BOHB Sampler") 

# Random
Tcr, poly_start, poly_degree, max_iter = bohb_opt_exp.minimizer;
opt_params_r = OptimizationParams(poly_start, poly_degree, max_iter);
# Parameters 
params = Parameters(grape_params, opt_params_r);
# Initial RF Pulse
B1ref = 1.0;
control_field = spline_RF(grape_params.N, Tcr, B1ref);
# Run Optimization
grape_output_ran = @time GrapeMR.grape(params, control_field, spins);
# Plots
plot_magnetization_control_field(grape_output_ran.control_field, grape_output_ran.isochromats)
plot_cost_values(grape_output_ran.cost_values, grape_output_ran.params.grape_params)

GrapeMR.run_cost_analysis(grape_output_ran.control_field, spins[end], 30.0, 20, grape_output_ran.params.grape_params.cost_function)

# BOHB
Tcb, poly_start, poly_degree, max_iter = bohb_opt_exp.minimizer;
opt_params_b = OptimizationParams(poly_start, poly_degree, 3500);
# Parameters 
params = Parameters(grape_params, opt_params_b);
# Initial RF Pulse
B1ref = 1.0;
control_field = spline_RF(grape_params.N, Tcb, B1ref);
# Run Optimization
grape_output_bohb = @time GrapeMR.grape(params, control_field, spins);
# Plots
plot_magnetization_control_field(grape_output_bohb.control_field, grape_output_bohb.isochromats)
plot_cost_values(grape_output_bohb.cost_values, grape_output_bohb.params.grape_params)


# Cost Function Analysis
GrapeMR.run_cost_analysis(grape_output_bohb.control_field, spins[end], 80.0, 20, grape_output_bohb.params.grape_params.cost_function)

# Save Optimization Data
save_grape_data(grape_output_ran; folder_path)
save_hyperopt_data(random_opt_exp; folder_path)

save_grape_data(grape_output_bohb; folder_path)
save_hyperopt_data(bohb_opt_exp; folder_path)

export_bruker(grape_output_ran)
export_bruker(grape_output_bohb)


iso_no_relax = dynamics.(grape_output_bohb.control_field, no_relax_spins) 
plot_transverse_time(iso_no_relax, grape_output_.control_field.t_control) 
