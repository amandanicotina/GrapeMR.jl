"""
    GrapeOutput{T, F}

Stores the output of GRAPE optimization.

# Fields
- `isochromats::Vector{Isochromat}`: Vector of isochromats containing the spin dynamics.
- `control_field::AbstractControlField`: Optimized control field (can be physical or normalized).
- `cost_values::Vector{Float64}`: Sequence of cost function values at each iteration.
- `params::Parameters{F}`: Struct containing GRAPE and optimization parameters.
"""
mutable struct GrapeOutput{T<:Real, CF<:AbstractControlField, F}
    isochromats::Vector{Isochromat}
    control_field::CF
    cost_values::Vector{Float64}
    params::Parameters{F}
end

function GrapeOutput(params::Parameters, control_field::AbstractControlField)
    T = eltype(control_field.B1x)
    config = params.opt_params.config
    max_iter = get_max_iter(config)
    return GrapeOutput{T, typeof(control_field), typeof(params.grape_params.cost_function)}(
        Isochromat[],
        deepcopy(control_field),
        zeros(T, max_iter),
        params
    )
end


