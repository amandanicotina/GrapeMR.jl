using GrapeMR
using TOML

tm = TOML.parsefile("src/default_config.toml")

offsets = collect(-tm["spins"]["offset"]:1:tm["spins"]["offset"])
spins = GrapeMR.Spin(
    tm["spins"]["M0"],
    [s["T1"] for s in tm["spins"]["intrinsics"]], 
    [s["T2"] for s in tm["spins"]["intrinsics"]],
    offsets, tm["spins"]["delta_B1"], 
    [s["target"] for s in tm["spins"]["intrinsics"]],
    [s["label"] for s in tm["spins"]["intrinsics"]]
)

# Grape Parameters
mask = tm["grape_parameters"]["fields2optimize"]
grape_params = GrapeParams(
    tm["grape_parameters"]["time_steps"], 
    Symbol(tm["grape_parameters"]["cost_function"]), 
    reshape(mask, 1, length(mask))  # we need a 1x3 Matrix{Bool} instead of a Vector{Bool}
)

# Optimization Parameters
opt_params   = OptimizationParams(
    tm["optimization_parameters"]["poly_start"], 
    tm["optimization_parameters"]["poly_degree"], 
    Int(ceil(tm["optimization_parameters"]["max_iter"]))
)

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
control_field = spline_RF(grape_params.N, tm["control_field"]["control_time"], tm["control_field"]["B1ref"]) 

# Run Optimization
grape_output = @time grape(params, control_field, spins); 

# Save data
if tm["save_files"]["enabled"]
    experiment_folder = save_grape_data(grape_output; folder_path= tm["save_files"]["folder_path"])
    if tm["save_files"]["export_bruker"]
        # Export data
        go = load_grape_data(experiment_folder)
        export_bruker(grape_output)
    end
end

if tm["plot"]
    # Plots
    plot_cost_values(grape_output.cost_values, grape_params)
    plot_magnetization_2D(grape_output.isochromats)
    plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
    plot_control_fields(grape_output.control_field; unit="Hz")
    plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control)
end

