using GrapeMR
using Test
using SafeTestsets

@testset "Check Bloch methods" begin include("test_bloch_methods.jl") end

@testset "Cost Function: Gradient vs Finite difference" begin include("test_cost_functions.jl") end

@testset "True gradient vs finite difference" begin include("test_cost_gradients.jl") end