using GrapeMR
using Test, SafeTestsets

@time begin
    @safetestset "Control Field Generation: construction, normalization, and shape types" begin
        include("control_fields/test_control_fields.jl")
    end

    @safetestset "Spin Data Type: struct definition and spin population constructor" begin
        include("bloch_methods/test_spin_generation.jl")
    end

    @safetestset "Bloch Simulation: normalized vs physical units (single spin)" begin
        include("bloch_methods/test_bloch_methods.jl")
    end

    @safetestset "Bloch Simulation: normalized vs physical units (multiple spins)" begin
        include("bloch_methods/test_bloch_methods_multi_spins.jl")
    end

    @safetestset "Bloch Simulation: Verification of Magnetization Trajectories" begin
        include("bloch_methods/test_bloch_methods_dynamics.jl")
    end    

    @safetestset "Cost Function Validation: Correctness and Gradient Consistency" begin
        include("cost_and_gradients/test_cost_functions_correctness.jl")
    end  

    @safetestset "Cost Function Gradient: Finite Differences vs Analytical" begin
        include("cost_and_gradients/test_cost_functions_correctness.jl")
    end    

    @safetestset "Gradient Computation: forward, adjoint, and field update steps" begin
        include("cost_and_gradients/test_gradients.jl")
    end

    # TODO 
    # @safetestset "Gradient Adjoint Check: Finite Differences vs Analytical" begin
    #     include("cost_and_gradients/test_gradient_correctness.jl")
    # end

    # @safetestset "grape() Convergence" begin include("test_grape.jl") end
    
    # @safetestset "Hyperparameter Optimization" begin include("optimization/test_optimize.jl") end
    
    # @safetestset "Tutorial" begin include("../docs/src/tutorial.jl") end
end
