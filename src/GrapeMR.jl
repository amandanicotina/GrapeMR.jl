module GrapeMR

using ArgParse
using BlochSim
using ColorSchemes
using CSV
using CubicSplines
using DataFrames
using Dates
using Distributed
using ForwardDiff
using Hyperopt
using JLD2
using LinearAlgebra
using Logging
using NumericalIntegration
using Optim
using ParameterSchedulers
using Plots
using PrettyPrint
using ProgressMeter
using Random
using StaticArrays
using TOML

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = SA[0 0 0 0; 0 0 0 0; 0 0 0 -1; 0 0 1 0]
const Iy = SA[0 0 0 0; 0 0 0 1; 0 0 0 0; 0 -1 0 0]

const γ_unit = 2π * γ_¹H     # [rad/s/T] — needed for internal Bloch matrices

# ----------- #
  # Includes #
# ----------- #

# Data types
include("data_types/ControlField.jl")
include("data_types/Spins.jl")
include("data_types/PulseShape.jl")
include("data_types/Optimizer.jl")
include("data_types/Parameters.jl")
include("data_types/GrapeOutput.jl")

# Pulse generators
include("rf_pulses/generators.jl")

# bSSFP module
include("bSSFP/data_types.jl")
include("bSSFP/steady_state.jl")
include("bSSFP/plots.jl")

# Analysis
include("analysis/rf_analysis.jl")
include("analysis/cost_analysis.jl")
include("analysis/magnetization_analysis.jl")

# Propagation and cost
include("propagation&cost/bloch_methods.jl")
include("propagation&cost/cost_functions.jl")

# Optimization
include("optimization/grape.jl")
include("optimization/gradients.jl")
include("optimization/hyperparameter_opt.jl")

# Utilities
include("utilities/save_data.jl")
include("utilities/export_bruker.jl")

# Plots
include("plots/plots_hyperparameters.jl")
include("plots/plots_control_field.jl")
include("plots/plots_magnetization.jl")


# ----------- #
  # Exports #
# ----------- #

# Constants
export γ_¹H, Ix, Iy

# Data types
export ControlField, NormalizedControlField
export Spins, Spin, Magnetization, Isochromat, generate_spins
export PulseShape, Spline, Hard, Sinc, Gaussian, BSSFP, pulse_shape
export AbstractOptimizer, GradientDescent, BFGS, ManualGradientDescent
export AbstractOptimizerConfig, GradientDescentConfig, BFGSConfig, ManualGradientDescentConfig
export OptimizationParams, GrapeParams, Parameters
export GrapeOutput


# RF pulse generation
export generate_control_field, normalize_control_field, denormalize_control_field

# GRAPE
export grape, grape!, dynamics, grape_gd_optim!
export backward_propagation!, forward_propagation!, bloch_matrix
export update!, gradient!, gradient

# Cost Functions
export euclidean_norm, spin_target, saturation_contrast

# File I/O
export save_grape_data, save_hyperopt_data, load_grape_data, load_hyperopt_data
export export_bruker

# Plotting
export plot_cost_values, plot_magnetization_control_field
export plot_control_fields, plot_control_fields_phase_shift
export plot_transverse_magnetization, plot_magnetization_2D, plot_magnetization_3D
export plot_magnetization_time, plot_transverse_time, plot_longitudinal_time
export initialize_plot, color_palette, get_target_properties

# Hyperparameter optimization
export random_hyperopt, bohb_hyperopt, hband_hyperopt
export plot_hyperopt_history, plot_cost_grape_runs, plot_hyperopt_contour, plot_cost_hyperparam
export plot_evaluations, plot_bohb

# ------------------ #
 # CLI Entrypoint #
# ------------------ #
function julia_main()::Cint
    s = ArgParseSettings()

    @add_arg_table! s begin
        "--config"
        help = "Path to the TOML configuration file."
        default = "src/default_config.toml"
    end

    parsed_args = parse_args(ARGS, s)
    run_grape_optimization(parsed_args["config"])
    return 0
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    GrapeMR.julia_main()
end
