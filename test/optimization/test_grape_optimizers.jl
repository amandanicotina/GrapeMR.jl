using GrapeMR
using Plots
using BenchmarkTools

# --- 1. Define comparison setup ---
# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = 0.0 #-15:1:15
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
## Gradient Descent
gd_config = GradientDescentConfig(0.75, 1, 500)
opt_params_gd = OptimizationParams(GradientDescent(), gd_config)

## Manual Gradient Descent
mgd_config = ManualGradientDescentConfig(0.75, 1, 500)
opt_params_mgd = OptimizationParams(ManualGradientDescent(), mgd_config)

## BFGS
bfgs_config = BFGSConfig(500)
opt_params_bfgs = OptimizationParams(BFGS(), bfgs_config)

# Combined parameter struct
params_gd = Parameters(grape_params, opt_params_gd)
params_mgd = Parameters(grape_params, opt_params_mgd)
params_bfgs = Parameters(grape_params, opt_params_bfgs)

# --- 2. Initial control field and deepcopy for each optimizer ---
control_field  = generate_control_field(:spline; t_c=0.5, B1ref=1.0)
cf_manual = deepcopy(control_field)
cf_gd     = deepcopy(control_field)
cf_bfgs   = deepcopy(control_field)

# --- 3. Initialize GrapeOutput structs ---
out_manual = GrapeOutput(params_gd, cf_manual)
out_gd     = GrapeOutput(params_mgd, cf_gd)
out_bfgs   = GrapeOutput(params_bfgs, cf_bfgs)

# --- 4. Run each optimizer with timing ---
println("Running manual gradient descent:")
@time grape!(out_manual, params_mgd, cf_manual, spins, ManualGradientDescent())

println("Running Optim GradientDescent():")
@time grape!(out_gd, params_gd, cf_gd, spins, GradientDescent())

println("Running Optim BFGS():")
@time grape!(out_bfgs, params_bfgs, cf_bfgs, spins, BFGS())

# --- 5. Compare cost function values ---
println("\nFinal cost values:")
@show out_manual.cost_values[end]
@show out_gd.cost_values[end]
@show out_bfgs.cost_values[end]

# --- 6. Plot cost evolution ---
p_cost = plot(
    [out_manual.cost_values, out_gd.cost_values, out_bfgs.cost_values],
    labels = ["Manual GD" "GD" "BFGS"],
    xlabel = "Iteration",
    ylabel = "Cost",
    lw = 2,
    title = "Cost Function Evolution"
)
display(p_cost)

# --- 7. Plot control field comparison ---
p_b1x = plot(cf_manual.B1x[1, :], label = "Manual GD", lw=1.5, xlabel = "Time step", ylabel = "Amplitude", title = "Optimized B1x")
plot!(p_b1x, cf_gd.B1x[1, :], label = "GD", lw=1.5)
plot!(p_b1x, cf_bfgs.B1x[1, :], label = "BFGS", lw=1.5)
display(p_b1x)

p_b1y = plot(cf_manual.B1y[1, :], label = "Manual GD", lw=1.5, xlabel = "Time step", ylabel = "Amplitude", title = "Optimized B1y")
plot!(p_b1y, cf_gd.B1y[1, :], label = "GD", lw=1.5)
plot!(p_b1y, cf_bfgs.B1y[1, :], label = "BFGS", lw=1.5)
display(p_b1y)

# # --- 7. Benchmark runtime (optional) ---
# println("\nBenchmarking (use @btime only if confident in precompilation):")
# @btime grape!($out_manual, $params_mgd, $cf_manual, $spins, ManualGradientDescent())
# @btime grape!($out_gd, $params_gd, $cf_gd, $spins, GradientDescent())
# @btime grape!($out_bfgs, $params_bfgs, $cf_bfgs, $spins, BFGS())


# @code_warntype grape(params_mgd, control_field, spins, ManualGradientDescent())
# @code_warntype grape(params_gd, control_field, spins, GradientDescent())
# @code_warntype grape(params_bfgs, control_field, spins, BFGS())

# @code_warntype grape!(out_manual, params_mgd, control_field, spins, ManualGradientDescent())
# @code_warntype grape!(out_gd, params_gd, control_field, spins, GradientDescent())
# @code_warntype grape!(out_bfgs, params_bfgs, control_field, spins, BFGS())
