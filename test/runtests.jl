using GrapeMR
using Test
using SafeTestsets

@testset "Check Bloch methods" begin include("test_bloch_methods.jl") end

@testset "Check true gradient vs finite difference" begin include("test_cost_gradients.jl") end