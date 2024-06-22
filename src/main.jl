using Plots
using GrapeMR
using JLD
using Hyperopt

# Grape Parameters 
grape_params = GrapeParams(1500, :target_one_spin, [true true false])

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

M0 = [0.0, 0.0, 1.0] 
T1 = [0.3]
T2 = [0.08]
label  = ["T1=800ms"]
target = ["max"]    
B0 = 50.0
offset = collect(0:1:B0) 
ΔB1 = [1.0]
spins = GrapeMR.Spin(M0, T1, T2, [0.0], ΔB1, target, label)

# Optimization Parameters
# bohb_max_iter = range(1000, stop=10000, step=500)
# bohb_Tc = LinRange(0.05, 1.0, 20)
# bohb = hyperoptimization(spins, grape_params, bohb_Tc, bohb_max_iter)
Tc, poly_start, poly_degree, max_iter = 1.0, 0.01, 1.0, 100 #   bohb.minimizer # 
opt_params   = OptimizationParams(poly_start, poly_degree, max_iter)

# RF
B1ref = 1.0
B1x = spline_RF(grape_params.N, Tc)'
B1y = spline_RF(grape_params.N, Tc)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Optimization
# grape_output = @time grape(opt_params, grape_params, control_field, spins)
rand = @time random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(100, step=100, stop=1000))
bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 20), range(100, step=100, stop=1000))
# bohb = @hyperopt for i = 50, sampler = Hyperband(R=50, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
#         Tc = Tc,
#         poly_start  = poly_start,
#         poly_degree = poly_degree,
#         max_iter    = max_iter;
   
#     if state !== nothing
#         Tc, poly_start, poly_degree, max_iter = state
#     end
#     if Tc >= 0.0 && poly_start >= 0.0 && poly_degree >= 1.0 && max_iter >= 1.0
#         print("\n", i, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t", max_iter, "\t")

#         # RFs
#         B1ref = 1.0
#         B1x = spline_RF(grape_params.N, Tc)'
#         B1y = spline_RF(grape_params.N, Tc)'
#         Bz  = zeros(1,grape_params.N)
#         control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

#         # Optimize
#         opt_params   = OptimizationParams(poly_start, poly_degree, max_iter)
#         res = grape(opt_params, grape_params, control_field, spins)

#         cost = res.cost_values[end]
#         cost, [Tc, poly_start, poly_degree, max_iter]
#     else
#         1000.0, [Tc, poly_start, poly_degree, max_iter]
#     end
# end


# Plots
# plot_magnetization_2D(grape_output.isochromats) 
# plot_control_fields(grape_output.control_field) 
# plot_cost_values(grape_output.cost_values, grape_params)
# plot_magnetization_time(grape_output.isochromats[2], Tc)

# pTrans = plot()
# isos = grape_output.isochromats
# for iso ∈ isos
#     m = iso.magnetization.dynamics
#     Mx = m[2,:]
#     My = m[3,:]
#     Mt = Mx + im*My
#    #plot!(pTrans, Mt, color = 2, label = false)
#     scatter!(pTrans, [Mt[end]], color = 2, label = false)
# end

# display(pTrans)
# # plot_Mtrans_offset_ss(grape_output.isochromats)
# # plot_magnetization_targetB0(grape_output.isochromats)
# # plot_magnetization_target_3D(grape_output.isochromats[5])

# # plot(bohb)  ==> issue > bohb.results[perm]
# # 1-element Vector{Any}: 
# # (0.7649095880816756, (0.6499999999999999, 0.01, 3))
# # This is the tuple and scatter cant't plot it like this

# # Save data
# # folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
# # GrapeMR.save_grape_data(grape_output; folder_path)

# ##### Phase testing 
# # ux = grape_output.control_field.B1x
# # uy = grape_output.control_field.B1y
# # ψ  = π
# # Ux = vec(ux*cos(ψ) .+ uy*sin(ψ)) 
# # Uy = vec(uy*cos(ψ) .- ux*sin(ψ))

# # tc = grape_output.control_field.t_control
# # t = range(0.0, tc, length = length(ux))
# # p_Bx = plot(t, Ux, linewidth=2, label=false, ylabel="B1x", title="Control Fields", titlefontsize=15)
# # p_By = plot(t, Uy, linewidth=2, label=false, ylabel="B1y", xlabel="t [s]")
# # # xticks_values = [-π, -π/2, 0, π/2, π]
# # # xticks_labels = ["-π", "-π/2", "0", "π/2", "π"]
# # # p_By = plot!(p_By, xticks=(xticks_values, xticks_labels))

# # p = plot(p_Bx, p_By, layout = (2,1))

# # plot(vec(Ux), vec(Uy))
# # using JLD
# # folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/"
# # folder_name = "2024-06-13/"
# # file_name   = "ss_target_1000ms_100Hz_5Hz.jl"
# # @save joinpath(folder_path, folder_name, file_name) grape_output opt_params

# # @save "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/2024-06-13/ss_target_1000ms_100Hz_5Hz.jl" grape_output opt_params




