"""
    OptimizationParams

Defines parameters for polynomial-based decay for ϵ in the GRAPE algorithm.

# Fields
- `poly_start::Float64`: Initial polynomial coefficient.
- `poly_degree::Int`: Degree of the polynomial.
- `max_iter::Int`: Maximum number of optimization iterations.
"""
struct OptimizationParams
    poly_start::Float64
    poly_degree::Int
    max_iter::Int
end

"""
    GrapeParams{F}

Encapsulates key parameters for the GRAPE optimization process.

# Fields
- `N::Int64`: Number of time steps in the control sequence.
- `cost_function::F`: Function used to compute the optimization cost.
- `fields_opt::Dict{String, Bool}`: Dictionary indicating which fields to optimize (e.g., `B1x`, `B1y`).
"""
struct GrapeParams{F}
    N::Int64
    cost_function::F
    fields_opt::Dict{String, Bool}
end

"""
    Parameters{F}

Combines GRAPE and optimization parameters for the full optimization process.

# Fields
- `grape_params::GrapeParams{F}`: GRAPE-specific optimization parameters.
- `opt_params::OptimizationParams`: General optimization parameters.
"""
struct Parameters{F}
    grape_params::GrapeParams{F}
    opt_params::OptimizationParams
end