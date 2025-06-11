using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = -15:1:15
T1 = [0.6, 0.3]
T2 = [0.1, 0.05]
label = ["C1", "C2"]
target = ["max", "min"]
spins = generate_spins(M0, T1, T2, offsets, ΔB1, target, label)

# GRAPE-specific parameters
grape_params = GrapeParams(
    GrapeMR.saturation_contrast,
    Dict("B1x" => true, "B1y" => true, "Bz" => false),
)

# Optimization parameters
gd_config = GradientDescentConfig(0.75, 1, 500)
opt_params_gd = OptimizationParams(GradientDescent(), gd_config)

# Manual Gradient Descent
mgd_config = ManualGradientDescentConfig(0.75, 1, 2000)
opt_params_mgd = OptimizationParams(ManualGradientDescent(), mgd_config)

# BFGS-specific config
bfgs_config = BFGSConfig(500)
opt_params_bfgs = OptimizationParams(BFGS(), bfgs_config)

# Combined parameter struct
params_gd = Parameters(grape_params, opt_params_gd)
params_mgd = Parameters(grape_params, opt_params_mgd)
params_bfgs = Parameters(grape_params, opt_params_bfgs)

# Initial control field
B1ref = 1.0
t_c = 0.5
control_field  = generate_control_field(:spline; t_c=t_c, B1ref=B1ref)

# Run Optimization
grape_output_gd = @time grape(params_gd, control_field, spins);
grape_output_mgd = @time grape(params_mgd, control_field, spins);
grape_output_bfgs = @time grape(params_bfgs, control_field, spins);

# Plots
plot_magnetization_control_field(grape_output_gd.control_field, grape_output_gd.isochromats)
plot_magnetization_control_field(grape_output_mgd.control_field, grape_output_mgd.isochromats)
plot_magnetization_control_field(grape_output_bfgs.control_field, grape_output_bfgs.isochromats)

plot_cost_values(grape_output_gd.cost_values, grape_output_gd.params.grape_params)
plot_cost_values(grape_output_mgd.cost_values, grape_output_mgd.params.grape_params)
plot_cost_values(grape_output_bfgs.cost_values, grape_output_bfgs.params.grape_params)


# output_gd = GrapeOutput(params_gd, control_field)
# output_mgd = GrapeOutput(params_mgd, control_field)
# output_bfgs = GrapeOutput(params_bfgs, control_field)

# grape!(output_gd, params_gd, control_field, spins, GradientDescent())
# grape!(output_mgd, params_mgd, control_field, spins, ManualGradientDescent())
# grape!(output_bfgs, params_bfgs, control_field, spins, BFGS())

# using Plots
# plot!(grape_output.isochromats[3].magnetization.dynamics')

# plot_hyperopt_history(random_opt; title = "Random Sampler") 
# plot_hyperopt_history(bohb_opt; title = "BOHB Sampler") 
# plot_cost_grape_runs(bohb_opt; plotlog = true, title = "BOHB Sampler")
# plot_cost_hyperparam(bohb_opt; grid_plot = true, title = "BOHB Sampler")

# Cost Function Analysis
# GrapeMR.run_cost_analysis(grape_output.control_field, spins[1], 80.0, 20, grape_output.params.grape_params.cost_function)

# Save Optimization Data
# folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
# save_grape_data(grape_output; folder_path)
# save_hyperopt_data(random_opt; folder_path)
