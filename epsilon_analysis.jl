
using GrapeMR, ColorSchemes, Plots


struct GrapeOutputEpsilon{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}, F}
    isochromats::Vector{Isochromat}
    control_field::ControlField{T, M1, Mz}
    cost_values::Vector{Float64}
    params::Parameters{F}
    epsilons::Vector{Float64}
end


function grape_fix_epsilon(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter + 1)
    cost_vals = zeros(eltype(cf.B1x), op.max_iter, 1)[:]
    epsilons  = zeros(Float64, op.max_iter, 1)[:]
    u1x = Matrix{eltype(cf.B1x)}(undef, 1, size(cf.B1x, 2))
    u1y = Matrix{eltype(cf.B1x)}(undef, 1, size(cf.B1x, 2))    
    grape_output = GrapeOutputEpsilon(Vector{Isochromat}(), deepcopy(cf), cost_vals, p, epsilons)    
    ∇x = zeros(eltype(cf.B1x), 1, gp.N)
    ∇y = zeros(eltype(cf.B1x), 1, gp.N)
    mag, adj = zeros(Float64, 4, gp.N + 1), zeros(Float64, 4, gp.N + 1)
    ϵ = 1e-3

    for (ϵ, i) ∈ zip(lr_scheduler, 1:op.max_iter)
    # for i ∈ 1:op.max_iter
        # ϵ   = max(ϵ, eps)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)
        grape_output.epsilons[i] = ϵ

        for spin ∈ spins
            # Forward Propagation 
            iso = dynamics(cf, spin)
            mag = iso.magnetization.dynamics
            cost_vars = gp.cost_function(iso)
            # Cost Variables
            cost, adj_ini = cost_vars
            grape_output.cost_values[i,1] += cost
            # Adjoint Propagation
            adj = backward_propagation(adj_ini, cf, iso)

            # Save Isochromats from the last iterations
            if i == op.max_iter
                push!(grape_output.isochromats, iso)
            end

            # Gradient
            if gp.fields_opt[1]
                ∇x .+= gradient(adj, mag, Ix)
            end
            if gp.fields_opt[2]
                ∇y .+= gradient(adj, mag, Iy)
            end
        end

        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x .= u1x
        cf.B1y .= u1y
    end

    grape_output.control_field.B1x .= u1x
    grape_output.control_field.B1y .= u1y

    # Print Infos
    final_cost = round(grape_output.cost_values[end], digits=3)
    # println("\n Final Cost Function Value = $final_cost \n")
    # RF_pulse_analysis(grape_output.control_field)

    return grape_output
end

function grape_decay_epsilon(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter + 1)
    cost_vals = zeros(eltype(cf.B1x), op.max_iter, 1)[:]
    epsilons  = zeros(Float64, op.max_iter, 1)[:]
    u1x = Matrix{eltype(cf.B1x)}(undef, 1, size(cf.B1x, 2))
    u1y = Matrix{eltype(cf.B1x)}(undef, 1, size(cf.B1x, 2))    
    grape_output = GrapeOutputEpsilon(Vector{Isochromat}(), deepcopy(cf), cost_vals, p, epsilons)    
    ∇x = zeros(eltype(cf.B1x), 1, gp.N)
    ∇y = zeros(eltype(cf.B1x), 1, gp.N)
    mag, adj = zeros(Float64, 4, gp.N + 1), zeros(Float64, 4, gp.N + 1)
    ϵ = 1e-3

    for (ϵ, i) ∈ zip(lr_scheduler, 1:op.max_iter)
    # for i ∈ 1:op.max_iter
        # ϵ   = max(ϵ, eps)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)
        grape_output.epsilons[i] = ϵ

        for spin ∈ spins
            # Forward Propagation 
            iso = dynamics(cf, spin)
            mag = iso.magnetization.dynamics
            cost_vars = gp.cost_function(iso)
            # Cost Variables
            cost, adj_ini = cost_vars
            grape_output.cost_values[i,1] += cost
            # Adjoint Propagation
            adj = backward_propagation(adj_ini, cf, iso)

            # Save Isochromats from the last iterations
            if i == op.max_iter
                push!(grape_output.isochromats, iso)
            end

            # Gradient
            if gp.fields_opt[1]
                ∇x .+= gradient(adj, mag, Ix)
            end
            if gp.fields_opt[2]
                ∇y .+= gradient(adj, mag, Iy)
            end
        end

        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x .= u1x
        cf.B1y .= u1y
    end

    grape_output.control_field.B1x .= u1x
    grape_output.control_field.B1y .= u1y

    # Print Infos
    final_cost = round(grape_output.cost_values[end], digits=3)
    # println("\n Final Cost Function Value = $final_cost \n")
    # RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


# Spin Parameters
M0     = [0.0, 0.0, 1.0]
ΔB1    = [1.0]
offset = -15:1:15
T1 = [0.6]# 1.0]
T2 = [0.3]#, 0.6]
label  = ["s1"]#, "test"]
target = ["Center of Bloch's Ball"]#, "min"]
spins  = Spin(M0, T1, T2, [0.0], ΔB1, target, label)

# Optimization Parameters
Tc = 1.0
poly_start  = 0.1
poly_degree = 1
max_iter    = 3000
opt_params  = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.euclidean_norm, [true true false])

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref) 

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins); 

# Plots
# color_palette(var::AbstractArray) = ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=length(var)))) 
# color_palette(var::Int) = ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=var))) 
# colors = color_palette(10) 

# cost_fix_ϵ = grape_output_fix.cost_values
# par_fix_ϵ  = grape_output_fix.params.grape_params
# cmin_fix_ϵ  = round(cost_fix_ϵ[end], digits=3)
# iter_fix_ϵ = range(0.0, stop=length(cost_fix_ϵ), length=length(cost_fix_ϵ))

# cost_poly_ϵ = grape_output_poly.cost_values
# par_poly_ϵ = grape_output_poly.params.grape_params
# cmin_poly_ϵ = round(cost_poly_ϵ[end], digits=4)
# iter_poly_ϵ = range(0.0, stop=length(cost_poly_ϵ), length=length(cost_poly_ϵ))

# pCost = plot(
#     # ylims = (0.0, 1.0),
#     xlabel="Iterations",
#     ylabel="Cost Value",
#     title="Cost Function Convergence",
#     guidefontsize=12,
#     legendfontsize=10,
#     tickfontsize=10,
#     titlefontsize=12,
#     framestyle=:box,
#     grid=false)
# plot!(pCost, iter_fix_ϵ, cost_fix_ϵ, label="ϵ = e-3", lw=2, color=colors[2])
# scatter!([iter_fix_ϵ[end]], [cost_fix_ϵ[end]], color=colors[2], markersize=4, label="Final Cost = $cmin_fix_ϵ")

# plot!(pCost, iter_poly_ϵ, cost_poly_ϵ, label="ϵ = polynomial decay", lw=2, color=colors[8])
# # scatter!([iter_poly_ϵ[end]], [cost_poly_ϵ[end]], color=colors[8], markersize=4, label="Final Cost = $cmin_poly_ϵ")

# fix_ϵ_vals =grape_output_fix.epsilons
# poly_ϵ_vals = grape_output_poly.epsilons

# pEpsilon = plot(
#     # ylims = (0.0, 1.0),
#     xlabel="Iterations",
#     ylabel="ϵ",
#     title="ϵ Value",
#     guidefontsize=12,
#     legendfontsize=10,
#     tickfontsize=10,
#     titlefontsize=12,
#     framestyle=:box,
#     grid=false)
# plot!(pEpsilon, iter, fix_ϵ_vals, label="Constant", lw=2, color=colors[2])
# plot!(pEpsilon, iter, poly_ϵ_vals, label="Polynomial Decay", lw=2, color=colors[8])


