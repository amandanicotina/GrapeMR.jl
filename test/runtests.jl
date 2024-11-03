using GrapeMR
using Test, SafeTestsets

@time begin
    @time @safetestset "Cost Function: Gradient vs Finite Differences" begin include("test_cost_functions.jl") end

    @time @safetestset "grape() Convergence" begin include("test_grape.jl") end
    
    # @time @safetestset "Hyperparameter Optimization" begin include("optimization/test_optimize.jl") end
   
    # TODO(anicotina): Fix tests with more time. Normalisation issue apparently.
    # @time @safetestset "Bloch Methods: BlochSim vs GrapeMR" begin include("test_bloch_methods.jl") end
    # @time @safetestset "Gradient vs Finite Differences" begin include("test_gradients.jl") end

end
