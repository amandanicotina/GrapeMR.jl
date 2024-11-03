using GrapeMR
using Test, SafeTestsets

@time begin
    @safetestset "Cost Function: Gradient vs Finite Differences" begin include("test_cost_functions.jl") end

    @safetestset "grape() Convergence" begin include("test_grape.jl") end

    @safetestset "Bloch Methods: BlochSim vs GrapeMR" begin include("test_bloch_methods.jl") end
    
    @safetestset "Hyperparameter Optimization" begin include("optimization/test_optimize.jl") end
    
    # @time @safetestset "Gradient vs Finite Differences" begin include("test_gradients.jl") end

    @safetestset "Tutorial" begin include("../docs/src/tutorial.jl") end
end
