module GrapeMR

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
const Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];

using CSV
using Plots
using Dates
using Match
using BlochSim
using DataFrames
using CubicSplines
using LinearAlgebra
using ParameterSchedulers
using ExponentialUtilities

include("data_types/ControlField.jl")
include("data_types/OptimizationParams.jl")
include("data_types/Spins.jl")

include("propagation&cost/bloch_methods.jl")
include("propagation&cost/steady_state.jl")
include("propagation&cost/cost_functions.jl")
include("propagation&cost/cost_gradients.jl")

include("set_optimization/optimize.jl")
include("set_optimization/finite_difference.jl")

include("utilities/save_data.jl")
include("utilities/rf_shapes.jl")

include("plots.jl")

export γ_¹H, Ix, Iy
export ControlField, ControlFieldbSSFP
export OptimizationParams
export Spins, Spin, SteadyState, Magnetization, Isochromat

export forward_propagation, backward_propagation
export steady_state, steady_state_matrix, steady_state_geometric
export cost_function
export cost_function_gradient

export grape, gradient, update!
export finite_difference_cost, finite_difference_field

export save_grape_data
export initial_field_spline

export plot_cost_values, plot_control_fields
export plot_magnetization_time, plot_magnetization_Mz_Mxy
export plot_magnetization_target, plot_magnetization_target_3d

end
