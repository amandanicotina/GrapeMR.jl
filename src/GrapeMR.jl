module GrapeMR

using ArgParse
using Plots
using CSV
using TOML
using JLD2
using Dates
using Distributed
using Logging
using ColorSchemes
using Hyperopt
using BlochSim
using DataFrames
using CubicSplines
using LinearAlgebra
using NumericalIntegration
using ParameterSchedulers
using PrettyPrint
using StaticArrays
using ForwardDiff

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = SA[0 0 0 0; 0 0 0 0; 0 0 0 -1; 0 0 1 0]
const Iy = SA[0 0 0 0; 0 0 0 1; 0 0 0 0; 0 -1 0 0]

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

include("plots/plots_hyperparameters.jl")
include("plots/plots_control_field.jl")
include("plots/plots_magnetization.jl")

export γ_¹H, Ix, Iy

# Data types
export ControlField
export OptimizationParams, GrapeParams, Parameters, GrapeOutput
export Spins, Spin, Magnetization, Isochromat, generate_spins

# Grape 
export gaussian_RF, spline_RF, sinc_RF, bSSFP_RF, hard_RF
export grape, dynamics, run_grape_optimization
export backward_propagation, backward_propagation!, forward_propagation
export finite_difference_cost, finite_difference_field, gradient

# Save/load/export Files
export save_grape_data, save_hyperopt_data, load_grape_data, load_hyperopt_data
export export_bruker

# Plots
export plot_cost_values, plot_magnetization_control_field
export plot_control_fields, plot_control_fields_phase_shift
export plot_transverse_magnetization, plot_magnetization_2D, plot_magnetization_3D
export plot_magnetization_time, plot_transverse_time, plot_longitudinal_time
export initialize_plot, color_palette, get_target_properties

# Hyperopt
export random_hyperopt, bohb_hyperopt, hband_hyperopt
export plot_hyperopt_history, plot_cost_grape_runs, plot_hyperopt_contour, plot_cost_hyperparam
export plot_evaluations, plot_bohb



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

end

# Calling the main function when the script is executed
if abspath(PROGRAM_FILE) == @__FILE__
    GrapeMR.julia_main()
end
