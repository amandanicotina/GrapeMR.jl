struct OptimizationParams
    poly_start::Float64
    poly_degree::Int
    max_iter::Int
end

struct GrapeParams{F}
    N::Int64
    cost_function::F
    fields_opt::Array{Bool, 2}
end

struct Parameters{F}
    grape_params::GrapeParams{F}
    opt_params::OptimizationParams
end

"""
List of cost function symbols
    :euclidean_norm      
    :target_one_spin    
    :steady_state      
    :steady_state_opt    
    :saturation_contrast 
    :saturation_contrast_square
"""
