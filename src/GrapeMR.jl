module GrapeMR

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
const Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];

using CSV
using JLD2
using Plots
using Dates
using Distributed
using Logging
using Match
using Hyperopt
using Symbolics
# using BlochSim
using DataFrames
using CubicSplines
using LinearAlgebra
using NumericalIntegration
using ParameterSchedulers
# using Wandb

include("data_types/ControlField.jl")
include("data_types/Parameters.jl")
include("data_types/Spins.jl")

include("analysis/rf_analysis.jl")
include("analysis/magnetization_analysis.jl")

include("propagation&cost/bloch_methods.jl")
include("propagation&cost/steady_state.jl")
include("propagation&cost/cost_functions.jl")

include("optimization/optimize.jl")
include("optimization/finite_difference.jl")

include("utilities/save_data.jl")
include("utilities/export_bruker.jl")
include("utilities/rf_shapes.jl")

include("plots/plots_bohb.jl")
include("plots/plots_control_field.jl")
include("plots/plots_cost_functions.jl")
include("plots/plots_magnetization.jl")
include("plots/plots_time_dynamics.jl")

export γ_¹H, Ix, Iy

export complex_signal, amplitudes_and_phases, bruker_normalized_amplitudes_and_phases
export integral_factor, fast_fourier_transform, average_pulse_power

export ControlField, ControlFieldbSSFP, ControlFieldNormalized
export OptimizationParams, GrapeParams, Parameters, GrapeOutput
export Spins, Spin, SpinNormalized, SteadyState, SpinRange, Magnetization, Isochromat

export forward_propagation, backward_propagation, bloch_matrix, normalization, inverse_normalization
export test_forward_propagation, test_backward_propagation
export steady_state, steady_state_matrix, steady_state_geometric, steady_state_geometric_Mz
export cost_function
export cost_function_gradient

export grape, par_grape, tpar_grape, norm_grape, random_sample, hyperoptimization
export gradient, update!
export finite_difference_cost, finite_difference_field, finite_difference_field_symmetric

export save_grape_data, load_grape_data, file_name_string
export export_bruker
export spline_RF, sinc_RF, bSSFP_RF, hard_RF

export plot_cost_values, plot_cost_offset
export plot_control_fields, plot_control_fields_phase, plot_control_fields_BxBy
export plot_magnetization_time, plot_magnetization_2D, plot_magnetization_3D
export plot_magnetization_target, plot_magnetization_target_3D
export plot_magnetization_targetB0, plot_Mtrans_offset_ss
export plot_transverse_magnetization, plot_magnetization_control_field
export plot_evaluations, plot_bohb

end