using Test
using GrapeMR
using SafeTestsets

@time begin
    @time @safetestset "Cost Function: Gradient vs Finite difference" begin 
        include("test_cost_functions.jl") 
    end
    @time @safetestset "Gradient vs Finite Differences" begin 
        include("test_gradients.jl") 
    end
end
