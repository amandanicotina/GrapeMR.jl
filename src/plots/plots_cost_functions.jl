
"""

Cost Function plots

"""

function plot_cost_values(cost::Vector{Float64}, gp::GrapeParams)
    p = plot(cost, label = string(gp.cost_function), lw = 2,
    ylims = (0.0, 1.0),
    xlabel = "Iterations",
    ylabel = "Cost Value",
    title  = "Cost Function Convergence",
    titlefontsize = 12,
    )

    return p
end

function plot_cost_offset(isos::Vector{Isochromat}, cost::Symbol)
    # Cost function values for all isochromats
    ans = cost_function.(isos, cost)
    c = [x[1] for x in ans]


    # Offset frequencies
    # ν_ini = isos[1].spin.B0inho
    # ν_end = isos[end].spin.B0inho
    # ν_len = Int(ceil(length(isos)/2))
    # ν = range(ν_ini, stop=ν_end, length=ν_len)

    # p = plot(xlabel = "Offset [Hz]",
    #     ylabel = "Cost Value",
    #     title  = "Cost Function Offset profile",
    #     titlefontsize = 12)
    #     plot!(p, ν, c[1:ν_len], label =  "min", lw = 2)
    #     plot!(p, ν, c[ν_len:end], label = "max", lw = 2)

    ν_ini = isos[1].spin.B0inho
    ν_end = isos[end].spin.B0inho
    ν_len = ceil(length(isos))
    ν = range(ν_ini, stop=ν_end, length=ν_len)
    p = plot(xlabel = "Offset [Hz]",
        ylabel = "Cost Value",
        title  = "Cost Function Offset profile",
        titlefontsize = 12)
        plot!(p, ν, c, label =  "min", lw = 2)
    return p
end

function plot_cost_offset(spins::Vector{GrapeMR.Spin}, cost::Symbol)
    # Calculate dynamics for new offset range with optimized field
    isos = Vector{Isochromat}()
    for spin ∈ spins
        mag = forward_propagation(grape_output.control_field, spin)
        dyn = GrapeMR.Magnetization(mag)
        iso = Isochromat(dyn, spin)
        push!(isos, iso)
    end
    # Cost function values for all isochromats
    ans = cost_function.(isos, cost)
    c = [x[1] for x in ans]

    # Offset frequencies
    ν_ini = isos[1].spin.B0inho
    ν_end = isos[end].spin.B0inho
    ν_len = ceil(length(isos))
    ν = range(ν_ini, stop=ν_end, length=ν_len)

    p = plot(xlabel = "Offset [Hz]",
        ylabel = "Cost Value",
        title  = "Cost Function Offset profile",
        titlefontsize = 12)
        plot!(p, ν, c, label =  "min", lw = 2)
        # plot!(p, ν, c[ν_len+1:end], label = "max", lw = 2)

    return p
end

# using GrapeMR, JLD2, Plots
# folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/"
# folder_name = "2024-06-10/"
# file_name   = "min_600ms_25Hz.jl"
# (go, op, gp) = JLD.@load joinpath(folder_path, folder_name, file_name) grape_output opt_params 

# M0 = [0.0, 0.0, 1.0]
# labeloff, targetoff  = ["T1 = 300ms"], ["min"]
# T1 = grape_output.isochromats[1].spin.T1
# T2 = grape_output.isochromats[1].spin.T2

# offset = 100  # [Hz]
# offset_vals  = collect(-offset:1:offset)
# B1_vals = collect(range(0.7, stop=1.3, length=5))

# spins_offset_B0 = GrapeMR.Spin(M0, T1, T2, offset_vals, [1.0], targetoff, labeloff)
# spins_offset_B1 = GrapeMR.Spin(M0, T1, T2, [0.0], B1_vals, targetoff, labeloff)
# spins_offset_B0_B1 = GrapeMR.Spin(M0, T1, T2, offset_vals, B1_vals, targetoff, labeloff)

# get_cost = :target_one_spin#grape_params.cost_function
# control_field = grape_output.control_field

# isos = Isochromat[]
# costs = [] 

# for spin in spins_offset_B0
#     mag = forward_propagation(grape_output.control_field, spin)
#     dyn = GrapeMR.Magnetization(mag)
#     iso = Isochromat(dyn, spin)
#     push!(isos, iso)
#     push!(costs, cost_function(iso, get_cost))
# end
# p = plot(offset_vals, costs, label=false, lw=2, color=:viridis,
#     ylims = [-0.05, 1.0],
#     xlabel = "Offset [Hz]",
#     ylabel = "Cost Value",
#     title  = "Cost Function Offset profile",
#     framestyle=:box,
#     titlefontsize = 12
# )

# cost_matrix = zeros(length(offset_vals), length(B1_vals))
# spins_matrix = [GrapeMR.Spin(M0, T1, T2, [b0], [b1], targetoff, labeloff) for b0 in offset_vals, b1 in B1_vals]
# for i in 1:length(offset_vals)
#     for j in 1:length(B1_vals)
#         spin = spins_matrix[i, j][]
#         mag = forward_propagation(control_field, spin)
#         dyn = GrapeMR.Magnetization(mag)
#         iso = Isochromat(dyn, spin)
#         cost_matrix[i, j] = cost_function(iso, get_cost)
#     end
# end

# c_m = cost_matrix/maximum(cost_matrix)
# h1 = contourf(B1_vals, offset_vals, c_m, color=:viridis, #framestyle=:box,
#      xlabel="B1 Inhomogeneity [%]", ylabel="Offset [Hz]", title="Cost Function Map", colorbar_title="Cost Value")