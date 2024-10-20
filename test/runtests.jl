using Test
using GrapeMR
using SafeTestsets


# TODO automatic differention test for calculating the cost funcs derivatives
@time begin
    # @time @safetestset "Cost Function: Gradient vs Finite difference" begin include("test_cost_functions.jl") end
    # @time @safetestset "Bloch Methods: BlochSim vs GrapeMR" begin include("test_bloch_methods.jl") end
    # @time @safetestset "Gradient vs Finite Differences" begin include("test_gradients.jl") end
    @safetestset "Optimization" begin include("optimization/test_optimize.jl") end
end
