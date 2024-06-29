using Plots
using GrapeMR
using JLD
using Hyperopt


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

# M0 = [0.0, 0.0, 1.0] 
# T1 = [0.3, 1.0]
# T2 = [0.08, 0.6]
# label  = ["T1=300ms", "T1=1000ms"]  
# target = ["max", "min"]
# B0 = 30.0
# offset = collect(-B0/2:1:B0/2) 
# ΔB1 = [1.0]
# spins = GrapeMR.Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
# bohb_max_iter = range(1000, stop=10000, step=500)
# bohb_Tc = LinRange(0.05, 1.0, 20)
# bohb = hyperoptimization(spins, grape_params, bohb_Tc, bohb_max_iter)
# Tc, poly_start, poly_degree, max_iter = 0.5, 0.1, 2, 2000 #   bohb.minimizer #
# opt_params   = OptimizationParams(poly_start, poly_degree, max_iter)

# Grape Parameters 
# grape_params = GrapeParams(2000, :target_phase_encoding, [true true false])

# RF
# B1ref = 1.0
# B1x = spline_RF(grape_params.N, Tc)'
# B1y = spline_RF(grape_params.N, Tc)'
# Bz  = zeros(1, grape_params.N)
# control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Optimization
# grape_output = @time grape(opt_params, grape_params, control_field, spins)
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 100), 5000, i=100)
# random_hopt = @time random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=100)

# Plots
# plot_magnetization_2D(grape_output.isochromats) 
# plot_control_fields(grape_output.control_field) 
# plot_cost_values(grape_output.cost_values, grape_params)
# plot_magnetization_time(grape_output.isochromats[2], Tc)

# pTrans = plot()#xlims=[-0.5, 0.5], ylims=[0.5, 1.0])
# isos = grape_output.isochromats
# for iso ∈ isos
#     m = iso.magnetization.dynamics
#     Mx = m[2,:]
#     My = m[3,:]
#     Mt = Mx + im*My
#    plot!(pTrans, Mt, color = 2, label = false)
#     scatter!(pTrans, [Mt[end]], color = 2, label = false)
# end

# display(pTrans)


# @save "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/2024-06-22/1spin_My_T1300ms_50Hz_phase.jl" grape_output opt_params grape_params
# plot_Mtrans_offset_ss(grape_output.isochromats)
# plot_magnetization_targetB0(grape_output.isochromats)
# plot_magnetization_target_3D(grape_output.isochromats[5])

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



M0 = [0.0, 0.0, 1.0]
T1 = [0.3, 1.0]
T2 = [0.08, 0.6]
label  = ["T1=300ms", "T1=1000ms"]  
target = ["max", "min"]
B0 = 30.0
offset = collect(-B0/2:1:B0/2) 
ΔB1 = [1.0]
spins = GrapeMR.Spin(M0, T1, T2, [0.0], ΔB1, target, label)

# Optimization Parameters
# bohb_max_iter = range(1000, stop=10000, step=500)
# bohb_Tc = LinRange(0.05, 1.0, 20)
# bohb = hyperoptimization(spins, grape_params, bohb_Tc, bohb_max_iter)
# Tc, poly_start, poly_degree, max_iter = 0.5, 0.1, 2, 2000 #   bohb.minimizer #
# opt_params   = OptimizationParams(poly_start, poly_degree, max_iter)

# Grape Parameters 
grape_params = GrapeParams(1000, :saturation_contrast, [true true false])

# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.0, 100), 5000, i=100)
rand_hopt = random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=20)


function bohb_params(bohb)
    cost  = []
    t_c   = []
    start = []
    deg   = []
    max   = []
    for i ∈ eachindex(bohb.results)
        push!(cost, bohb.results[i][1])
        push!(t_c, bohb.history[i][1])
        push!(start, bohb.history[i][2])
        push!(deg, bohb.history[i][3])
        push!(max, bohb.history[i][4])
    end
    function get_eps(start, deg, max_iter)
        return start / (1 - (max_iter/2 - 1) / max_iter)^deg
    end
    eps = [get_eps(start, deg, it) for (start, deg, it) in zip(start, deg, max)]
    t_min = bohb.minimizer[1]
    st_min = bohb.minimizer[2]
    d_min = bohb.minimizer[3]
    m_min = bohb.minimizer[4]
    c_min = round(bohb.minimum[1], digits=4)
    order = collect(1:length(t_c))
    m_iter = 1000:100:5000
    p_cost =  scatter(t_c, max, zcolor=cost, 
         markerstrokecolor = :auto, label = false,
         xlabel = "Control time [s]", ylabel = "Iter", colorbar_title="Cost Value",
         title = "Optimization",
         color = :viridis)
         scatter!([t_min], [m_min], label = "Minimum = $c_min", 
         marker = :star5, markersize = 8, color = :red)
    p_order =  scatter(cost, m_iter, zcolor=order, 
         markerstrokecolor = :auto, label = false, 
         xlabel = "Cost Function", ylabel = "Resources - log scale", colorbar_title="Iterations",
         title = "Hyperparameter Tuning - Random Sampler", 
         color = :viridis)
         scatter!([c_min], [m_min], label = "Minimum = $c_min", 
         marker = :star5, markersize = 8, color = :red)

    return p_cost, p_order
end

p1, p2 = bohb_params(rand_hopt)
p1
p2

rand_hopt

using JLD2 
JLD2.@save "/Users/daviddodel/git-repos/GrapeMR.jl/rand_hopt.jld2" rand_hopt
JLD2.@save "/Users/daviddodel/Documents/amanda/rand_hopt.jl" rand_hopt
JLD.@save "/Users/daviddodel/Documents/amanda/rand_hopt.jld" rand_hopt
JLD.@save "/Users/daviddodel/Documents/amanda/rand_hopt.jl" rand_hopt



