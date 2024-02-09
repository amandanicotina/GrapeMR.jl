module GrapeMR

using Plots
using Zygote
using LinearAlgebra
using ExponentialUtilities

include("get_params.jl")
include("plots.jl")

include("propagation/bloch_methods.jl")

include("set_optimization/optimize.jl")
include("set_optimization/cost_functions.jl")
include("set_optimization/cost_gradients.jl")
include("set_optimization/get_gradient.jl")

export γ_¹H 
export InitialControlFields, Spins, OptimizationParams, Magnetization
export plot_magnetization, plot_cost_values

export forward_propagation, backward_propagation

export grape_optimize, finite_difference_cost, finite_difference_field

export cost_functions, cost_gradients
export normalization, gradient_controls, update_control_field

end
