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
include("set_optimization/finite_difference.jl")

include("set_optimization/cost_functions.jl")
include("set_optimization/cost_gradients.jl")


#include("plots.jl")

export γ_¹H, Ix, Iy
export ControlField
export OptimizationParams
export Spin, Magnetization, Isochromat

export normalization, forward_propagation, backward_propagation

export grape_optimize, grape
export finite_difference_cost, finite_difference_field

export euclidean_norm, target_one_spin, saturation_contrast
export grad_euclidean_norm, grad_target_one_spin, grad_saturation_contrast
export gradient, update

#export plot_magnetization, plot_cost_values, plot_control_fields, plot_magnetization_target

end
