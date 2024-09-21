# using Distributed
# const num_processes::Int = 6;
# addprocs(num_processes)
# println(nprocs()) 
# println(workers())

# @everywhere begin
using GrapeMR

# Spin System
## Gycloproteins
# T1_gycons = collect(range(650, stop = 1100, length=2).*1e-3)
# T2_gycons = collect(range(10, stop = 100, length=2).*1e-3)
# label_gycons  = fill("Gycloproteins", length(T1_gycons)) 
# target_gycons = fill("max", length(T1_gycons))
## Water
# T1_water = collect(range(1150, stop = 3000, length=2).*1e-3)
# T2_water = collect(range(1500, stop = 3000, length=2).*1e-3)
# label_water  = fill("Water", length(T1_water)) 
# target_water = fill("min", length(T1_water))
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = 50.0
target = ["[0.0, 1.0, 0.0]"]
label  = ["s1"]
T1 = [1.0]
T2 = [0.6]
offset = collect(-B0:5:B0) 
# spins = GrapeMR.Spin(M0, [T1_gycons; T1_water], [T2_gycons; T2_water], B0, ΔB1, [target_gycons; target_water], [label_gycons; label_water])
spins = GrapeMR.Spin(M0, T1, T2, [-1.0, 0.0, 1.0], ΔB1, target, label)

# Optimization Parameters
Tc = 0.8
poly_start = 0.1
poly_degree = 1.0
max_iter = 750.0
opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(2000, :euclidean_norm, [true true false])

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output = @time grape(params, control_field, spins); 
#norm_grape_output = @time norm_grape(params, control_field, spins); 

# Save data
folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
experiment_folder = save_grape_data(grape_output; folder_path)
go = load_grape_data(experiment_folder)

# export_bruker(res_grape)

# Plots
grape_output.cost_values[end]
plot_cost_values(grape_output.cost_values/maximum(grape_output.cost_values), grape_params)
plot_magnetization_2D(grape_output.isochromats)
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control)


# using Plots
# mag_analysis = false
# if mag_analysis == true

# using DataFrames, Plots
# function topspin_dynamics_comparison(file_path::String, iso::GrapeMR.Isochromat; s::String="Longitudinal")
#     # Call the function and store the result
#     rf_time = 1.0
#     data = read_nmrsim_data(file_path)
#     if s == "Transverse"
#         p = plot_tmagnetization_simul_vs_topspin(data, iso, rf_time)
#     else 
#         p = plot_magnetization_simul_vs_topspin(data, iso, rf_time) 
#     end

#     return p
# end

# function read_nmrsim_data(file_path::String)
#     # Empty arrays to hold the data
#     Mx = Float64[]
#     My = Float64[]
#     Mz = Float64[]

#     # Open file
#     open(file_path, "r") do file
#         for line in eachline(file)
#             if startswith(line, ";")
#                 continue
#             end
#             # Split the line into components by whitespace and convert to Float64
#             values = split(strip(line), r"\s+")
#             if length(values) == 3
#                 push!(Mx, parse(Float64, values[1]))
#                 push!(My, parse(Float64, values[2]))
#                 push!(Mz, parse(Float64, values[3]))
#             end
#         end
#     end

#     # Return the data as a DataFrame
#     return DataFrame(Mx=Mx, My=My, Mz=Mz)
# end


# function plot_magnetization_simul_vs_topspin(df::DataFrame, iso::GrapeMR.Isochromat, t::Float64)
#         # TopSpin nmrsim
#         Mx_ts = df.Mx
#         My_ts = df.My
#         Mz_ts = df.Mz
    
#         # Grape
#         m  = iso.magnetization.dynamics
#         s = iso.spin
#         time = range(0.0, t, length = length(m[1,:]))

#         pMx = plot(time, m[2,:], label = false, lw = 2, color = cgrad(:viridis)[128],
#             ylabel = "Mx",
#             title  = "Magnetization Dynamics - Sample = $(s.label), Target = $(s.target)",
#             titlefontsize = 12)
#             scatter!(time, Mx_ts, label = false, markersize = 0.8, color = "black")
    
#         pMy = plot(time, -m[3,:], label = false, lw = 2, color = cgrad(:viridis)[128],
#             ylabel = "My")
#             scatter!(time, My_ts, label = false, markersize = 0.8, color = "black")

    
#         pMz =   plot(time, m[4,:], label = "GrapeMR", lw = 2, color = cgrad(:viridis)[128], 
#                 ylabel = "Mz",
#                 xlabel = "t [s]",
#                 legend = :bottomleft)
#                 scatter!(time, Mz_ts, label = "TopSpin", markersize = 0.8, color = "black")
#         pMag = plot(pMx, pMy, pMz, layout = (3,1))
    
#         return pMag
# end

# function plot_tmagnetization_simul_vs_topspin(df::DataFrame, iso::GrapeMR.Isochromat, t::Float64)
#     # TopSpin nmrsim
#     Mz_ts = df.Mz
#     Mt_ts = sqrt.(df.Mx.^2 + df.My.^2)

#     # Grape
#     m  = iso.magnetization.dynamics
#     s = iso.spin
#     time = range(0.0, t, length = length(m[1,:]))
#     mtrans = sqrt.(m[2,:].^2 + m[3,:].^2)

#     pTrans = plot(time, mtrans, label = false, lw = 2, color = cgrad(:viridis)[128],
#                 ylabel = "Transverse",                
#                 xlabel = "t [s]",
#                 legend = :bottomleft)
#                 scatter!(time, Mt_ts, label = false, markersize = 1, color = "black")
#     pMz =   plot(time, m[4,:], label = "GrapeMR", lw = 2, color = cgrad(:viridis)[128], 
#                 ylabel = "Longitudinal",
#                 title  = "Magnetization Dynamics - Sample = $(s.label)", #, Target = $(s.target)",
#                 titlefontsize = 12)

#                 scatter!(time, Mz_ts, label = "TopSpin", markersize = 1, color = "black")
#     pMag = plot(pMz, pTrans, layout = (2,1))

#     return pMag
# end

# spins = GrapeMR.Spin([0.0, 0.0, 1.0], [1e8], [1e8], [0.0], [1.0], ["-"], ["-"])
# mag = forward_propagation(grape_output.control_field, spins[1])
# dyn = GrapeMR.Magnetization(mag)
# iso = Isochromat(dyn, spins[1])
# filepath="/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/Metabolomics/TopSpin/target_one_spin_26-08.txt"
# topspin_dynamics_comparison(filepath, iso)

# end



# using Profile
# @profview grape(opt_params, grape_params, control_field, spins)

# using Plots

# function plot_grad_x(g)
#     px=plot() 
#     for i in 1:length(g[1])
#         grad_x = vec(g[1][i])
#         plot!(px,grad_x, label = false, title = "∇x")
#     end
#     return px
# end
# function plot_grad_y(g)
#     py=plot() 
#     for i in 1:length(g[2])
#         grad_y = vec(g[2][i])
#         plot!(py, grad_y, label = false, title = "∇y")
#     end
#     return py
# end

# px = plot_grad_x(g)
# py = plot_grad_y(g)

# display(px)  
# display(py) 

# # Hyperparameter searct
# # const wandb_project::String = "GrapeMR"
# # logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, n ame = nothing)
# # bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.5, 100), 5000, i=100)
# # rand_hopt = random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=20)
# # close(logger)
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.1, 1.5, 100), 3000, i=10)