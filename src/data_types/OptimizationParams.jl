
struct OptimizationParams
    N::Int64
    max_iter::Int64
    cost_function::Function # □ use @enum to get all the cost functions
    fields_opt::Array{Bool, 2}
end

