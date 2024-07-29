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
T1 = [0.8, 0.08]
T2 = [0.2, 0.04]
label  = ["T1=$(T1[1]*1e3)ms", "T1=$(T1[2]*1e3)ms"]  
target = ["max", "min"]
B0 = 30.0
offset = collect(-B0/2:3:B0/2) 
ΔB1 = [1.0]
spins = Spin(M0, T1, T2, [0.0], ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, :saturation_contrast_Mx, [true true false])

# Optimization Parameters
#bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 100), 5000, i=10)
Tc, poly_start, poly_degree, max_iter = 0.836, 0.1, 2, 100 # bohb.minimizer # 
opt_params   = OptimizationParams(poly_start, poly_degree,  Int(ceil(max_iter)))

# RF
B1ref = 1.0
B1x = spline_RF(grape_params.N, Tc)'
B1y = spline_RF(grape_params.N, Tc)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Optimization
grape_output = @time grape(opt_params, grape_params, control_field, spins)

# JLD2.@save joinpath(folder_path, folder_name, file_name) grape_output opt_params grape_params

# # Plots
# grape_output.cost_values[end]
# plot_cost_values(grape_output.cost_values, grape_params)
# plot_magnetization_2D(grape_output.isochromats)
# plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)

