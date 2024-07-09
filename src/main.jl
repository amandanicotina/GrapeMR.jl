using Plots
using GrapeMR
using JLD
using Hyperopt
using Profile

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
T1 = [1.3, 0.73]
T2 = [2.0, 0.04]
label  = ["T1=$(T1[1]*1e3)ms", "T1=$(T1[2]*1e3)ms"]  
target = ["min", "max"]
B0 = 30.0
offset = collect(-B0/2:5:B0/2) 
ΔB1 = [1.0]
spins = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 100), 5000, i=10)
Tc, poly_start, poly_degree, max_iter = 0.836, 0.1, 1, 2000 #  bohb.minimizer # 0
opt_params   = OptimizationParams(poly_start, poly_degree,  Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, :saturation_contrast_Mx, [true true false])

# RF
B1ref = 1.0
B1x = spline_RF(grape_params.N, Tc)'
B1y = spline_RF(grape_params.N, Tc)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Optimization
grape_output_euromar = @time grape(opt_params, grape_params, control_field, spins);
grape_output_test = @time par_grape(opt_params, grape_params, control_field, spins)

# grape_output = @time test_grape(opt_params, grape_params, control_field, spins);

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
