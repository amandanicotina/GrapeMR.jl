module GrapeMR

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
const Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];


using Plots
using Zygote
using SumTypes
using LinearAlgebra
using ExponentialUtilities

include("data_types/ControlField.jl")
include("data_types/OptimizationParams.jl")
include("data_types/Spins.jl")

include("propagation/bloch_methods.jl")

include("set_optimization/optimize.jl")

include("set_optimization/cost_functions.jl")
include("set_optimization/cost_gradients.jl")
include("set_optimization/get_gradient.jl")

include("plots.jl")

export γ_¹H, Ix, Iy
export ControlFields
export OptimizationParams
export Spins, Magnetization

export normalization, forward_propagation, backward_propagation

export grape_optimize, finite_difference_cost, finite_difference_field

export cost_functions, cost_gradients
export get_gradient, update_control_field

export plot_magnetization, plot_cost_values, plot_control_fields

end
