using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = 0.0
T1 = [1.0, 0.25] #[1/31.3436]
T2 = [0.08, 0.04] #[1/37.6471]
label = ["C1", "C2"]
target = ["min", "max"]
spins = generate_spins(M0, T1, T2, offsets, ΔB1, target, label)

# GRAPE-specific parameters
N = 1000
Tc = 0.5
B1ref = 5.0

grape_params = GrapeParams(
    N,
    GrapeMR.euclidean_norm,
    Dict("B1x" => true, "B1y" => true, "Bz" => false),
)

# Optimization parameters
# random_opt = @time random_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), range(2500, 5000, step = 500)) 
# bohb_opt = @time bohb_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 3000)
# hband_opt = @time hband_hyperopt(spins, grape_params, LinRange(0.1, 0.5, 10), 10)
poly_start, poly_degree, max_iter = 0.75, 1, 5000
opt_params = OptimizationParams(poly_start, poly_degree, max_iter)

# Combined parameter struct
params = Parameters(grape_params, opt_params)

# Initial control field generation (refactored using trait-based dispatch)
control_field  = generate_control_field(:spline; t_c=Tc, B1ref=B1ref)

# Run Optimization
grape_output = @time grape(params, control_field, spins)


# Plots
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)


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
