"""
    GrapeOutput{T, M1, Mz, F}

Stores the output of GRAPE optimization.

# Fields
- `isochromats::Vector{Isochromat}`: Vector of isochromats containing the spin dynamics.
- `control_field::ControlField{T, M1, Mz}`: Optimized control field after GRAPE optimization.
- `cost_values::Vector{Float64}`: Sequence of cost function values at each iteration.
- `params::Parameters{F}`: Struct containing GRAPE parameters and optimization parameters.
"""
mutable struct GrapeOutput{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}, F}
    isochromats::Vector{Isochromat}
    control_field::ControlField{T, M1, Mz}
    cost_values::Vector{Float64}
    params::Parameters{F}
end