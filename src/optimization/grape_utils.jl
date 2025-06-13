get_max_iter(config::GradientDescentConfig) = config.max_iter
get_max_iter(config::ManualGradientDescentConfig) = config.max_iter
get_max_iter(config::BFGSConfig) = config.max_iter

"""
    GrapeContext{T, CF, S}

Context container for shared data used during GRAPE optimization.

# Fields
- `control_field::CF`: The control field being optimized (e.g., B1x/B1y matrices).
- `spins::Vector{S}`: The ensemble of spin systems.
- `fields_opt::Dict{String, Bool}`: Dictionary specifying which control fields to optimize (`"B1x"` and/or `"B1y"`).
- `grape_params::GrapeParams`: Parameters for GRAPE, including cost function and optimization flags.
- `magnetization::Matrix{T}`: Reusable array for storing magnetization trajectories.
- `adjoint::Matrix{T}`: Reusable array for storing adjoint variables during backpropagation.
"""
mutable struct GrapeContext{T, CF, SV, F}
    control_field::CF
    spins::SV
    fields_opt::Dict{String, Bool}
    grape_params::GrapeParams{F}
    magnetization::Matrix{T}
    adjoint::Matrix{T}
end

"""
    build_grape_context(control_field, spins, grape_params) -> GrapeContext

Initializes and returns a `GrapeContext` object with preallocated memory
for magnetization and adjoint vectors.

# Arguments
- `control_field`: The initial RF control field (e.g., of type `AbstractControlField`).
- `spins`: A vector of spin systems (e.g., `Vector{Spin}`).
- `grape_params`: The `GrapeParams` struct that includes the cost function and field optimization flags.

# Returns
- `GrapeContext`: A reusable context for GRAPE optimization iterations.
"""
function build_grape_context(control_field, spins, grape_params::GrapeParams{F}) where F
    n = size(control_field.B1x, 2)
    T = eltype(control_field.B1x)
    return GrapeContext{T, typeof(control_field), typeof(spins), F}(
        control_field,
        spins,
        grape_params.fields_opt,
        grape_params,
        zeros(T, 4, n + 1),
        zeros(T, 4, n + 1),
    )
end

"""
    cost_function_vec!(u_vec, ctx::GrapeContext) -> Float64

Evaluates the cost function for a flattened control vector `u_vec`, using the provided `ctx`.

# Arguments
- `u_vec`: Flattened vector containing the control field values (`B1x`, `B1y`).
- `ctx::GrapeContext`: Shared context with control field, spin systems, and cost function.

# Returns
- `Float64`: Total cost for the control field configuration.
"""
function cost_function_vec!(u_vec::Vector{Float64}, ctx::GrapeContext)
    n = size(ctx.control_field.B1x, 2)
    if ctx.fields_opt["B1x"]
        ctx.control_field.B1x .= reshape(view(u_vec, 1:n), 1, :)
    end
    if ctx.fields_opt["B1y"]
        ctx.control_field.B1y .= reshape(view(u_vec, n+1:2n), 1, :)
    end

    total_cost = 0.0
    for spin in ctx.spins
        forward_propagation!(ctx.magnetization, ctx.control_field, spin)
        iso = Isochromat(Magnetization(ctx.magnetization), spin)
        cost, _ = ctx.grape_params.cost_function(iso)
        total_cost += cost
    end
    return total_cost
end

"""
    gradient_function_vec!(G, u_vec, ctx::GrapeContext) -> AbstractVector

Computes the gradient of the cost function with respect to the control field vector `u_vec`.

# Arguments
- `G`: Preallocated gradient vector to store the result.
- `u_vec`: Flattened control field vector to evaluate at.
- `ctx::GrapeContext`: Context with control fields, spin systems, and optimization settings.

# Returns
- `G`: Modified in place with the gradient vector.
"""
function gradient_function_vec!(G, u_vec::Vector{Float64}, ctx::GrapeContext)
    n = size(ctx.control_field.B1x, 2)
    if ctx.fields_opt["B1x"]
        ctx.control_field.B1x .= reshape(view(u_vec, 1:n), 1, :)
    end
    if ctx.fields_opt["B1y"]
        ctx.control_field.B1y .= reshape(view(u_vec, n+1:2n), 1, :)
    end

    T = eltype(ctx.magnetization)
    grad_x = zeros(T, 1, n)
    grad_y = zeros(T, 1, n)    

    for spin in ctx.spins
        forward_propagation!(ctx.magnetization, ctx.control_field, spin)
        iso = Isochromat(Magnetization(ctx.magnetization), spin)
        _, adj_init = ctx.grape_params.cost_function(iso)
        backward_propagation!(ctx.adjoint, ctx.control_field, iso, adj_init)

        if ctx.fields_opt["B1x"]
            grad_x .+= gradient(ctx.adjoint, ctx.magnetization, Ix)
        end
        if ctx.fields_opt["B1y"]
            grad_y .+= gradient(ctx.adjoint, ctx.magnetization, Iy)
        end
    end

    G .= vcat(vec(grad_x), vec(grad_y))
    return G
end


function early_stopping(res::Optim.MultivariateOptimizationResults)
    @info "check for custom stopping"
    return Optim.converged(res)
end
