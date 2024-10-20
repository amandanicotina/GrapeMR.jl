module GrapeMR

const γ_¹H = 42.5774688e6 #[Hz/T] 
const Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
const Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];

using CSV
using JLD2
# using Plots
using Dates
using Distributed
using Logging
# using ColorSchemes
using Match
# using Hyperopt
using Symbolics
using BlochSim
using DataFrames
using CubicSplines
using LinearAlgebra
using NumericalIntegration
using ParameterSchedulers
# using Wandb

include("data_types/ControlField.jl")
include("data_types/Parameters.jl")
include("data_types/Spins.jl")

include("bSSFP/data_types.jl")
include("bSSFP/steady_state.jl")
# include("bSSFP/plots.jl")

include("analysis/rf_analysis.jl")
# include("analysis/cost_analysis.jl")
# include("analysis/magnetization_analysis.jl")

include("propagation&cost/bloch_methods.jl")
include("propagation&cost/cost_functions.jl")

include("optimization/optimize.jl")
include("optimization/opt_test_func.jl")
include("optimization/finite_difference.jl")

include("utilities/save_data.jl")
include("utilities/export_bruker.jl")
include("utilities/rf_shapes.jl")

# include("plots/plots_bohb.jl")
# include("plots/plots_control_field.jl")
# include("plots/plots_magnetization.jl")

export γ_¹H, Ix, Iy

# Data types
export ControlField, ControlFieldbSSFP
export OptimizationParams, GrapeParams, Parameters, GrapeOutput
export Spins, Spin, SpinRange, Magnetization, Isochromat

# Analysis
# export complex_signal, amplitudes_and_phases, bruker_normalized_amplitudes_and_phases
# export integral_factor, average_pulse_power # check exportinh these functions when the run_rf_analysis is ready

# export run_cost_analysis

# 
export SteadyState, SteadyStateData
# export calculate_steady_state, plot_ss_offset_profile, plot_ss_flip_angle_profile
export steady_state, steady_state_matrix, steady_state_geometric, steady_state_geometric_Mz  # check exportinh these functions when the run_rf_analysis is ready

# Grape 
export grape, random_sample, hyperoptimization, dynamics
export threads_grape, old_grape, no_threads_grape

# Save/load/export files
export save_grape_data, save_bohb_data, load_grape_data, load_bohb_data
export export_bruker
export spline_RF, sinc_RF, bSSFP_RF, hard_RF

# Plots
# export plot_cost_values, plot_magnetization_control_field
# export plot_control_fields, plot_control_fields_phase_shift
# export plot_magnetization_time, plot_magnetization_2D, plot_magnetization_3D
# export plot_magnetization_target, plot_magnetization_target_3D
# export plot_magnetization_targetB0, plot_Mtrans_offset_ss
# export plot_transverse_magnetization, plot_transverse_time
# export plot_evaluations, plot_bohb


function julia_main()::Cint
    println("$ARGS")

    M0 = [0.0, 0.0, 1.0]
    ΔB1 = [1.0]
    B0 = 5.0
    offsets = collect(-B0:1:B0)

    # Water
    T1_water = 0.5
    T2_water = 0.1
    label_water = "S1"
    target_water = "[0.0, 1.0, 0.0]"

    # Glycerol
    spins = GrapeMR.Spin(M0, [T1_water], [T2_water], offsets, ΔB1, [target_water], [label_water])

    # Grape Parameters 
    # grape_params = GrapeParams(1500, :spin_target, [true true false])
    grape_params = GrapeParams(1500, "spin_target")

    # Optimization Parameters
    # bohb = @time hyperoptimization(spins, grape_params, LinRange(0.01, 1.0, 15), 1500)
    # plot_bohb(bohb)
    # spline_bohb = @time hyperoptimization(spins, grape_params, LinRange(0.01, 1.0, 15), 1500)
    # Tc, poly_start, poly_degree, max_iter = spline_bohb.minimizer
    Tc = 0.5
    poly_start = 0.1
    poly_degree = 1.0
    max_iter = 10.0
    # max_iter = 1000.0
    opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

    # Parameters 
    params = Parameters(grape_params, opt_params)

    # Initial RF Pulse
    B1ref = 1.0
    control_field = spline_RF(grape_params.N, Tc, B1ref) 

    # Run Optimization
    grape_output = grape(params, control_field, spins);
    return 0 # if things finished successfully
end

end