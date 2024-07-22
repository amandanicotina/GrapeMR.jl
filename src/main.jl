using Plots
using GrapeMR

# using JLD
# using Hyperopt
# using Profile

# Spin parameters
# M0 = [0.0, 0.0, 1.0] 
# T1 = [1.83, 0.62]#, 2.430];
# T2 = [0.18, 0.1]#, 0.132];
# label  = ["Blastoderm", "Yolk"]#, "White"];
# target = ["max", "min"]#, "min"];
# B0 = 10.0        
# offset = collect(-B0:2:B0) 
# ΔB1 = [1.0]
# # α, Δϕ, TR, TE = 2π/9, π, 5e-3, 5e-3/2
# # spins = GrapeMR.SteadyState(M0, T1, T2, offset, ΔB1, target, label, α, Δϕ, TR, TE)
# spins = GrapeMR.Spin(M0, T1, T2, offset, ΔB1, target, label)

# t1water = 1.3
# t2water = 2.0
# t2gycon = 0.04
# t1gycon = 0.730

M0 = [0.0, 0.0, 1.0] 
T1 = [1.3]#, 0.73]
T2 = [2.0]#, 0.04]
label  = ["T1=$(T1[1]*1e3)ms"]#, "T1=$(T1[2]*1e3)ms"]  
target = ["min"]#, "max"]
B0 = 30.0
offset = collect(-B0/2:3:B0/2) 
ΔB1 = [1.0]
spins = Spin(M0, T1, T2, [0.0], ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, :euclidean_norm, [true true false])

# Optimization Parameters
#bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 100), 5000, i=10)
Tc, poly_start, poly_degree, max_iter = 0.836, 0.1, 1, 1000 # bohb.minimizer # 
opt_params   = OptimizationParams(poly_start, poly_degree,  Int(ceil(max_iter)))

# RF
B1ref = 1.0
B1x = spline_RF(grape_params.N, Tc)'
B1y = spline_RF(grape_params.N, Tc)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Optimization
grape_output_test1 = @time grape(opt_params, grape_params, control_field, spins)
# 139.937298 seconds (2.62 G allocations: 174.438 GiB, 5.13% gc time, 0.00% compilation time)
# 177.664103 seconds (2.62 G allocations: 174.526 GiB, 7.60% gc time)
grape_output_test2 = @time par_grape(opt_params, grape_params, control_field, spins);
# 146.689483 seconds (2.62 G allocations: 174.561 GiB, 6.65% gc time, 0.32% compilation time)
# 169.531075 seconds (2.62 G allocations: 174.477 GiB, 8.51% gc time)
grape_output_test3 = @time tpar_grape(opt_params, grape_params, control_field, spins);
# 457.104520 seconds (2.62 G allocations: 174.545 GiB, 2.49% gc time, 0.09% compilation time)
# 469.271114 seconds (2.62 G allocations: 174.476 GiB, 2.72% gc time)
grape_output_test4 = @time GrapeMR.optimized_grape(opt_params, grape_params, control_field, spins);
# 447.482420 seconds (2.62 G allocations: 174.909 GiB, 2.74% gc time, 0.08% compilation time)

# random_hopt = @time random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=2)

# using JLD
# folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/"
# folder_name = "2024-06-27/"
# file_name   = "euromar_max_T180ms_30Hz.jl"
# (go, op, gp) = JLD.@save joinpath(folder_path, folder_name, file_name) grape_output opt_params grape_params
# JLD.@save joinpath(folder_path, folder_name, file_name) grape_output opt_params grape_params

# Plots
plot_cost_values(grape_output_euromar.cost_values, grape_params)
plot_magnetization_2D(grape_output_euromar.isochromats)
plot_magnetization_control_field(grape_output_euromar.control_field, grape_output_euromar.isochromats)

# plot_magnetization_2D(grape_output.isochromats) 
# plot_transverse_magnetization(grape_output_euromar.isochromats) 
# plot_control_fields(grape_output.control_field) 
# plot_cost_values(grape_output.cost_values, grape_params)
# plot_magnetization_time(grape_output_euromar.isochromats[1], Tc)

###### Save data
# folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
# GrapeMR.save_grape_data(grape_output; folder_path)

###### Test new backward_propagation function
# iso = grape_output.isochromats[1]
# cost_grad = GrapeMR.cost_function_gradient(iso, grape_params.cost_function)
# @time adj0 = backward_propagation(control_field, iso, cost_grad)
# @time adj1 = test_backward_propagation(control_field, iso, cost_grad)

# plot(adj0[1,:])
# plot!(adj1[1,:])

# plot(adj0[2,:])
# plot!(adj1[3,:])

# plot(adj0[3,:])
# plot!(adj1[3,:])

# plot(adj0[4,:])
# plot!(adj1[4,:])
