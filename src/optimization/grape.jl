get_max_iter(config::GradientDescentConfig) = config.max_iter
get_max_iter(config::ManualGradientDescentConfig) = config.max_iter
get_max_iter(config::BFGSConfig) = config.max_iter


function cost_function_vec!(u_vec, control_field, spins, magnetization, cost_function, fields_opt)
    n = size(control_field.B1x, 2)
    if fields_opt["B1x"]
        control_field.B1x .= reshape(view(u_vec, 1:n), 1, :)
    end
    if fields_opt["B1y"]
        control_field.B1y .= reshape(view(u_vec, n+1:2n), 1, :)
    end

    total_cost = 0.0
    for spin in spins
        forward_propagation!(magnetization, control_field, spin)
        iso = Isochromat(Magnetization(magnetization), spin)
        cost, _ = cost_function(iso)
        total_cost += cost
    end
    return total_cost
end

function gradient_function_vec!(G, u_vec, control_field, spins, magnetization, adjoint, cost_function, fields_opt)
    n = size(control_field.B1x, 2)
    if fields_opt["B1x"]
        control_field.B1x .= reshape(view(u_vec, 1:n), 1, :)
    end
    if fields_opt["B1y"]
        control_field.B1y .= reshape(view(u_vec, n+1:2n), 1, :)
    end

    grad_x = zeros(Float64, 1, n)
    grad_y = zeros(Float64, 1, n)

    for spin in spins
        forward_propagation!(magnetization, control_field, spin)
        iso = Isochromat(Magnetization(magnetization), spin)
        _, adj_init = cost_function(iso)
        backward_propagation!(adjoint, control_field, iso, adj_init)

        if fields_opt["B1x"]
            grad_x .+= gradient(adjoint, magnetization, Ix)
        end
        if fields_opt["B1y"]
            grad_y .+= gradient(adjoint, magnetization, Iy)
        end
    end

    G .= vcat(vec(grad_x), vec(grad_y))
    return G
end

"""
    grape(params, control_field, spins, optimizer)

Run the GRAPE algorithm using the specified optimizer.
Returns a `GrapeOutput` with cost values and optimized control fields.
"""
function grape(params::Parameters{F}, control_field::AbstractControlField, spins::Vector{<:Spins}, optimizer::GradientDescent) where F
    output = GrapeOutput(params, control_field)
    return grape!(output, params, output.control_field, spins, optimizer)
end
function grape(params::Parameters{F}, control_field::AbstractControlField, spins::Vector{<:Spins}, optimizer::BFGS) where F
    output = GrapeOutput(params, control_field)
    return grape!(output, params, output.control_field, spins, optimizer)
end
function grape(params::Parameters{F}, control_field::AbstractControlField, spins::Vector{<:Spins}, optimizer::ManualGradientDescent) where F
    output = GrapeOutput(params, control_field)
    return grape!(output, params, output.control_field, spins, optimizer)
end

function grape(params::Parameters, control_field::AbstractControlField, spins::Vector{<:Spins})
    grape(params, control_field, spins, params.opt_params.optimizer)
end

#########################################
# GRAPE Implementation - Gradient Descent
#########################################

"""
    grape!(output, params, control_field, spins, opt::GradientDescent)

In-place GRAPE optimizer using Gradient Descent.

# Returns
- `output`: Updated with optimized isochromats, cost values, and control field.
"""
function grape!(output::GrapeOutput, params::Parameters, control_field::AbstractControlField, spins::Vector{<:Spins}, opt::ManualGradientDescent)
    opt_params_config, grape_params = params.opt_params.config, params.grape_params
    scheduler = Poly(start=opt_params_config.poly_start, degree=opt_params_config.poly_degree, max_iter=opt_params_config.max_iter + 1)

    n = size(control_field.B1x, 2)
    T = eltype(control_field.B1x)

    ∇x = zeros(T, 1, n)
    ∇y = zeros(T, 1, n)
    u1x = similar(control_field.B1x, 1, n)
    u1y = similar(control_field.B1y, 1, n)

    fill!(output.cost_values, 0.0)
    empty!(output.isochromats)

    @showprogress for (ϵ, i) in zip(scheduler, 1:opt_params_config.max_iter)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)

        for spin in spins
            magnetization = zeros(Float64, 4, n + 1)
            adjoint = zeros(Float64, 4, n + 1)
            forward_propagation!(magnetization, control_field, spin)
            dyn = Magnetization(magnetization)
            iso = Isochromat(dyn, spin)
        
            cost, adj_init = grape_params.cost_function(iso)
            output.cost_values[i] += cost
            backward_propagation!(adjoint, control_field, iso, adj_init)
        
            if i == length(output.cost_values)
                push!(output.isochromats, iso)
            end
        
            if grape_params.fields_opt["B1x"]
                ∇x .+= gradient(adjoint, magnetization, Ix)
            end
            if grape_params.fields_opt["B1y"]
                ∇y .+= gradient(adjoint, magnetization, Iy)
            end
        end

        u1x, u1y = update!(control_field, (∇x, ∇y), ϵ)
        control_field.B1x .= u1x
        control_field.B1y .= u1y
        c = output.cost_values[i]
        # @info "Cost Function Value = $(round(c, digits=3))"
    end

    @info "Final Cost Function Value = $(round(output.cost_values[end], digits=3))"
    @debug "Average cost over iterations: $(mean(output.cost_values))"

    return output
end

#########################################
# GRAPE Implementation - BFGS
#########################################

"""
    grape!(output, params, control_field, spins, opt::BFGS)

In-place GRAPE optimizer using BFGS.
# Returns
- `output`: Updated with optimized isochromats, cost values, and control field.
"""
function grape!(output::GrapeOutput, params::Parameters, control_field::AbstractControlField, spins::Vector{<:Spins}, opt::BFGS)
    opt_params_config, grape_params = params.opt_params.config, params.grape_params
    fields_opt = grape_params.fields_opt
    n = size(control_field.B1x, 2)

    output.cost_values .= NaN
    magnetization = zeros(Float64, 4, n + 1)
    adjoint = zeros(Float64, 4, n + 1)

    u0_parts = Vector{Vector{Float64}}()
    if fields_opt["B1x"]
        push!(u0_parts, vec(control_field.B1x))
    end
    if fields_opt["B1y"]
        push!(u0_parts, vec(control_field.B1y))
    end
    u0 = reduce(vcat, u0_parts)

    max_iter = get_max_iter(opt_params_config)
    pbar = Progress(max_iter, desc = "Running GRAPE $(opt) Optimizer")
    iteration_counter = Ref(0)

    function progress_callback(state)
        iteration_counter[] += 1
        i = iteration_counter[]
        if i ≤ length(output.cost_values)
            output.cost_values[i] = state.value
        end
        next!(pbar)
        @debug "Iter: $i | Cost: $(round(state.value, digits=5))"
        return false
    end

    options = Optim.Options(
        iterations = max_iter,
        show_trace = false,
        callback = progress_callback,
        show_every = 1,
    )

    result = Optim.optimize(
        u_vec -> cost_function_vec!(u_vec, control_field, spins, magnetization, grape_params.cost_function, fields_opt),
        (G, u_vec) -> gradient_function_vec!(G, u_vec, control_field, spins, magnetization, adjoint, grape_params.cost_function, fields_opt),
        u0,
        LBFGS(),
        options
    )

    u_opt = Optim.minimizer(result)

    if fields_opt["B1x"]
        control_field.B1x .= reshape(view(u_opt, 1:n), 1, :)
    end
    if fields_opt["B1y"]
        control_field.B1y .= reshape(view(u_opt, n+1:2n), 1, :)
    end

    final_cost = cost_function_vec!(u_opt, control_field, spins, magnetization, grape_params.cost_function, fields_opt)
    output.cost_values[end] = final_cost

    empty!(output.isochromats)
    for spin in spins
        forward_propagation!(magnetization, control_field, spin)
        iso = Isochromat(Magnetization(copy(magnetization)), spin)
        push!(output.isochromats, iso)
    end
    @info "Final Cost Function Value = $(round(output.cost_values[end], digits=3))"

    return output
end

function grape!(output::GrapeOutput, params::Parameters, control_field::AbstractControlField, spins::Vector{<:Spins}, opt::GradientDescent)
    opt_params_config, grape_params = params.opt_params.config, params.grape_params
    fields_opt = grape_params.fields_opt
    n = size(control_field.B1x, 2)

    output.cost_values .= NaN
    magnetization = zeros(Float64, 4, n + 1)
    adjoint = zeros(Float64, 4, n + 1)

    u0_parts = Vector{Vector{Float64}}()
    if fields_opt["B1x"]
        push!(u0_parts, vec(control_field.B1x))
    end
    if fields_opt["B1y"]
        push!(u0_parts, vec(control_field.B1y))
    end
    u0 = reduce(vcat, u0_parts)

    max_iter = get_max_iter(opt_params_config)
    pbar = Progress(max_iter, desc = "Running GRAPE $(opt) Optimizer")
    iteration_counter = Ref(0)

    function progress_callback(state)
        iteration_counter[] += 1
        i = iteration_counter[]
        if i ≤ length(output.cost_values)
            output.cost_values[i] = state.value
        end
        next!(pbar)
        @debug "Iter: $i | Cost: $(round(state.value, digits=5))"
        return false
    end

    options = Optim.Options(
        iterations = max_iter,
        show_trace = false,
        callback = progress_callback,
        show_every = 1,
    )

    result = Optim.optimize(
        u_vec -> cost_function_vec!(u_vec, control_field, spins, magnetization, grape_params.cost_function, fields_opt),
        (G, u_vec) -> gradient_function_vec!(G, u_vec, control_field, spins, magnetization, adjoint, grape_params.cost_function, fields_opt),
        u0,
        Optim.GradientDescent(),
        options;
        inplace = true,
        autodiff = :false
    )    

    u_opt = Optim.minimizer(result)

    if fields_opt["B1x"]
        control_field.B1x .= reshape(view(u_opt, 1:n), 1, :)
    end
    if fields_opt["B1y"]
        control_field.B1y .= reshape(view(u_opt, n+1:2n), 1, :)
    end

    final_cost = cost_function_vec!(u_opt, control_field, spins, magnetization, grape_params.cost_function, fields_opt)
    output.cost_values[end] = final_cost

    empty!(output.isochromats)
    for spin in spins
        forward_propagation!(magnetization, control_field, spin)
        iso = Isochromat(Magnetization(copy(magnetization)), spin)
        push!(output.isochromats, iso)
    end
    @info "Final Cost Function Value = $(round(output.cost_values[end], digits=3))"

    return output
end