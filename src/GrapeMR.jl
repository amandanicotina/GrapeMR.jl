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
using ColorSchemes
using Hyperopt
using Symbolics
using BlochSim
using DataFrames
using CubicSplines
using LinearAlgebra
using NumericalIntegration
using ParameterSchedulers

include("data_types/ControlField.jl")
include("data_types/Parameters.jl")
include("data_types/Spins.jl")

include("bSSFP/data_types.jl")
include("bSSFP/steady_state.jl")
include("bSSFP/plots.jl")

include("analysis/rf_analysis.jl")
include("analysis/cost_analysis.jl")
include("analysis/magnetization_analysis.jl")

include("propagation&cost/bloch_methods.jl")
include("propagation&cost/cost_functions.jl")

include("optimization/optimize.jl")
include("optimization/hyperparameter_opt.jl")
include("optimization/finite_difference.jl")

include("utilities/save_data.jl")
include("utilities/export_bruker.jl")
include("utilities/rf_shapes.jl")

include("plots/plots_bohb.jl")
include("plots/plots_control_field.jl")
include("plots/plots_magnetization.jl")

export γ_¹H, Ix, Iy

# Data types
export ControlField
export OptimizationParams, GrapeParams, Parameters, GrapeOutput
export Spins, Spin, SpinRange, Magnetization, Isochromat

# bSSFP
export SteadyState, SteadyStateData
export calculate_steady_state, plot_ss_offset_profile, plot_ss_flip_angle_profile
export steady_state, steady_state_matrix, steady_state_geometric, steady_state_geometric_Mz 

# Analysis
export complex_signal, amplitudes_and_phases, bruker_normalized_amplitudes_and_phases
export integral_factor, fast_fourier_transform, average_pulse_power # check export in these functions when the run_rf_analysis is ready
export run_cost_analysis 

# Grape 
export grape, dynamics, backward_propagation
export random_sampler, bohb_hyperopt, hyperband_hyperopt
export finite_difference_cost, finite_difference_field

# Save/load/export files
export save_output_data, load_grape_data, load_bohb_data
export export_bruker
export spline_RF, sinc_RF, bSSFP_RF, hard_RF

# Plots
export plot_cost_values, plot_magnetization_control_field
export plot_control_fields, plot_control_fields_phase_shift
export plot_magnetization_time, plot_magnetization_2D, plot_magnetization_3D
export plot_magnetization_target, plot_magnetization_target_3D
export plot_magnetization_targetB0, plot_Mtrans_offset_ss
export plot_transverse_magnetization, plot_transverse_time
export plot_evaluations, plot_bohb

end