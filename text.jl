t = Template(;
           user="amandanicotina",
           authors=["Amanda Nicotina", "David Dodel"],
           plugins=[
               License(name="MIT"),
               Git(),
               GitHubActions(),
           ],
       )



##########################################################################################


no_threads = [295.8, 295.9, 289.6, 297.7, 294.3]
multi_threads = [191.5 212.4 190.7 194.0 191.7]
speed_up = [1.54 1.39 1.52 1.53 1.54]

using DataFrames

# Create the DataFrame
data = DataFrame(
    TimeSingleThread = [295.8, 295.9, 289.6, 297.7, 294.3],
    TimeMultiThread = [191.5, 212.4, 190.7, 194.0, 191.7],
    SpeedUp = [1.54, 1.39, 1.52, 1.53, 1.54])

# Display the table
println(data)
using Plots

# Data
input_sizes = [1000, 5000, 10000, 25000, 50000]
time_single_thread = [295.8, 295.9, 289.6, 297.7, 294.3]
time_multi_thread = [191.5, 212.4, 190.7, 194.0, 191.7]
speed_up = [1.54, 1.39, 1.52, 1.53, 1.54]

# Plot Execution Time
p1 = plot(time_single_thread, label="Single Thread", marker=:o, ylabel="Execution Time (s)", title="Execution Time Comparison", legend=:topleft)
plot!(p1, time_multi_thread, label="Multi Thread", marker=:o)

# Plot Speedup
p2 = plot(speed_up, label="Speedup", marker=:o, ylabel="Speedup", title="Speedup Achieved", legend=:topleft)

# Display the plots
plot(p1, p2, layout=(1, 2), size=(1000, 400))






############################# MAGNETIZATIOpTrans = plot()
isos = grape_output3.isochromats
for iso ∈ isos
    m = iso.magnetization.dynamics
    Mt =  m[2,:] + im*m[3,:]
    Mx = real.(Mt)
    My = imag.(Mt)
    plot!(pTrans, Mx, My, color = "black", label = false, 
    title = "Transverse Magnetization",
    titlefontsize = 12,
        xlabel="Mxy", ylabel="Mz")
    scatter!(pTrans, [Mx[end]], [My[end]], color = "black", label = false, framestyle=:box)
end
pMag = plot_magnetization_2D(grape_output3.isochromats)
cf = grape_output3.control_field
time = range(0.0, cf.t_control, length = length(cf.B1x))*1e3
Bx = cf.B1x
By = cf.B1y
B1 = vec(Bx + im * By)

p_Bx = plot(time, abs.(B1), linewidth=2, grid = false, label=false, color=:viridis, 
        framestyle=:box, ylabel="|B1| [Hz]", title="Control Fields", titlefontsize=12)
p_By = plot(time, angle.(B1), linewidth=2, grid = false, label=false, color=:viridis, 
        framestyle=:box, ylabel="ϕ [rad]", xlabel="t [ms]", 
        yticks=([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "π/2", "π"]))

l = @layout [[a{0.7h}; b{0.3h}] c]
pMagCf = plot(p_Bx, p_By, pMag, layout = l)
display(pMagCf)
###################################################################################################### 





############################# B0 AND B1 ANALYSIS ##############################  
using GrapeMR, JLD2, Plots
folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/"
folder_name = "2024-06-27/"
file_name   = "euromar_max_T180ms_30Hz.jl"
(go, op, gp) = JLD2.@load joinpath(folder_path, folder_name, file_name) grape_output3 opt_params grape_params

M0 = [0.0, 0.0, 1.0]
T1 = [0.08, 1.0]
T2 = [0.03, 0.6]
labeloff  = ["T1=80ms", "T1=1000ms"]  
targetoff = ["max", "min"]
# T1 = [0.08]
# T2 = [0.03]
# labeloff  = [ "T1=1000ms"]  
# targetoff = ["min"]

offset = 50  # [Hz]
offset_vals  = collect(-offset:0.5:offset)
B1_vals = collect(range(0.7, stop=1.3, length=20))

spins_offset_B0 = GrapeMR.Spin(M0, T1, T2, offset_vals, [1.0], targetoff, labeloff)
spins_offset_B1 = GrapeMR.Spin(M0, T1, T2, [0.0], B1_vals, targetoff, labeloff)
spins_offset_B0_B1 = GrapeMR.Spin(M0, T1, T2, offset_vals, B1_vals, targetoff, labeloff)

get_cost = :saturation_contrast
control_field = grape_output2.control_field

isos = Isochromat[]
costs = [] 

    for spin in spins_offset_B0
        mag = forward_propagation(control_field, spin)
        dyn = GrapeMR.Magnetization(mag)
        iso = Isochromat(dyn, spin)
        push!(isos, iso)
        push!(costs, cost_function(iso, get_cost))
    end

    costsmul = costs*100
    norm_cost = costsmul/maximum(costsmul)
    # Offset frequencies
    ν_ini = isos[1].spin.B0inho
    ν_end = isos[end].spin.B0inho
    ν_len = Int(ceil(length(isos)/2))
    ν = range(ν_ini, stop=ν_end, length=ν_len)

    p = plot(label = false, lw=2,
        #ylims = [-0.05, 1.0],
        xlabel = "Offset [Hz]",
        ylabel = "Cost Value",
        title  = "Cost Function Offset profile",
        framestyle=:box,
        titlefontsize = 12
    )
        plot!(p, offset_vals, norm_cost[1:ν_len], label =  "min", lw = 2, color=:viridis)
        plot!(p, offset_vals, norm_cost[ν_len+1:end], label = "max", lw = 2, color = cgrad(:viridis)[128] )


cost_matrix = zeros(length(offset_vals), length(B1_vals))
spins_matrix = [GrapeMR.Spin(M0, [0.08], [0.03], [b0], [b1], targetoff, labeloff) for b0 in offset_vals, b1 in B1_vals]
for i in 1:length(offset_vals)
    for j in 1:length(B1_vals)
        spin = spins_matrix[i, j][]
        mag = forward_propagation(control_field, spin)
        dyn = GrapeMR.Magnetization(mag)
        iso = Isochromat(dyn, spin)
        cost_matrix[i, j] = cost_function(iso, get_cost)
    end
end

c_m = cost_matrix/maximum(cost_matrix)
h1 = contourf(B1_vals, offset_vals, c_m, color=:viridis, #framestyle=:box,
     xlabel="B1 Inhomogeneity [%]", ylabel="Offset [Hz]", title="Cost Function Map", colorbar_title="Cost Value")
######################################################################################################


################## BOHB ################## 
using JLD2, Plots, JLD
JLD2.@load "/Users/amandanicotina/Documents/Documents/Euromar2024/bohb.jld2" bohb 
JLD2.@load "/Users/amandanicotina/Downloads/random_sampler.jld2" random_sampler 

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
    # pos_cost   = (1 .+ cost)
    norm_cost  = cost/maximum(cost)
    norm_c_min = (1 + c_min)/maximum(norm_cost)

    p_cost = scatter(t_c, log.(max),  zcolor=(1 .+ cost), 
         markerstrokecolor = :auto, label = false,
         xlabel = "Control time [s]", ylabel = "Resources", colorbar_title="Cost Value",
         title = "Optimization",
         color = :viridis)
         scatter!([t_min],[log(m_min)],  label = "Min = $c_min", 
         legend=:bottomleft,
         marker = :star5, markersize = 8, color = :red)
    p_order =  scatter(cost, log.(max), zcolor=order, 
         markerstrokecolor = :auto, label = false, 
         xlabel = "Cost Function", ylabel = "Resources - log scale", colorbar_title="Iterations",
         title = "Hyperparameter Tuning - Random Sampler", 
         color = :viridis)
         scatter!([c_min], [log(m_min)], label = "Min = $c_min", 
         legend=:bottomleft,
         marker = :star5, markersize = 8, color = :red)

    return p_cost, p_order, norm_cost
end

(p1, p2, m) = bohb_params(bohb);
# display(p1)
display(p2)
(p3, p4, m) = bohb_params(random_sampler);
display(p3)
display(p4)
m













