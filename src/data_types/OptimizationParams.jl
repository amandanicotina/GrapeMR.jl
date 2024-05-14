
struct OptimizationParams
    N::Int64
    max_iter::Int64
    cost_function::Symbol # □ use @enum to get all the cost functions
    fields_opt::Array{Bool, 2}
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