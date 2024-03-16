module GrapeMR

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
const Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];


using Plots
using Dates
using Flux
using CubicSplines
using Serialization
using LinearAlgebra
using ParameterSchedulers
using ExponentialUtilities

include("data_types/ControlField.jl")
include("data_types/OptimizationParams.jl")
include("data_types/Spins.jl")

include("propagation&cost/bloch_methods.jl")
include("propagation&cost/cost_functions.jl")
include("propagation&cost/cost_gradients.jl")

include("set_optimization/optimize.jl")
include("set_optimization/finite_difference.jl")
include("set_optimization/initialization.jl")

include("utilities/save_data.jl")

include("plots.jl")

export γ_¹H, Ix, Iy
export ControlField
export OptimizationParams
export Spin, Magnetization, Isochromat

export forward_propagation, backward_propagation, steady_state
export cost_euclidean_norm, cost_target_one_spin, cost_saturation_contrast, saturation_contrast_square
export grad_euclidean_norm, grad_target_one_spin, grad_saturation_contrast, grad_saturation_contrast_square

export grape, gradient, update
export finite_difference_cost, finite_difference_field

export normalization, initial_field_spline

export save_grape_data

export plot_cost_values, plot_control_fields
export plot_magnetization_time, plot_magnetization_Mz_Mxy

end
