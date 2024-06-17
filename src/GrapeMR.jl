module GrapeMR

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
const Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];

using CSV
using JLD
using Plots
using Dates
using Match
using Hyperopt
using BlochSim
using DataFrames
using CubicSplines
using LinearAlgebra
using ParameterSchedulers

include("data_types/ControlField.jl")
include("data_types/Parameters.jl")
include("data_types/Spins.jl")

include("propagation&cost/bloch_methods.jl")
include("propagation&cost/steady_state.jl")
include("propagation&cost/cost_functions.jl")
include("propagation&cost/cost_gradients.jl")

include("optimization/optimize.jl")
include("optimization/finite_difference.jl")

include("utilities/save_data.jl")
include("utilities/rf_shapes.jl")

include("plots.jl")

export γ_¹H, Ix, Iy

export ControlField, ControlFieldbSSFP
export OptimizationParams, GrapeParams
export Spins, Spin, SteadyState, Magnetization, Isochromat

export forward_propagation, backward_propagation
export steady_state, steady_state_matrix, steady_state_geometric, steady_state_geometric_Mz
export cost_function
export cost_function_gradient

export grape, par_grape, hyperopt_grape
export gradient, update!
export finite_difference_cost, finite_difference_field

export save_grape_data
export spline_RF, sinc_RF, bSSFP_RF

export plot_cost_values, plot_cost_offset, plot_control_fields
export plot_magnetization_time, plot_magnetization_2D, plot_magnetization_3D
export plot_magnetization_target, plot_magnetization_target_3D
export plot_magnetization_targetB0, plot_Mtrans_offset_ss

end
