
struct OptimizationParams
    N::Int64
    cost_function::String # □ use @enum to get all the cost functions
    fields_opt::Array{Bool, 2}
end

