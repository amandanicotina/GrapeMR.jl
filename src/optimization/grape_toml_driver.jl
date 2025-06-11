

"""
    run_grape_optimization(config_path::String)

Runs the GRAPE optimization process based on a TOML configuration file. Initializes spins, sets up parameters, executes optimization, and optionally saves and plots the results.

# Arguments
- `config_path::String`: Path to the TOML configuration file containing parameters for spins, optimization, and control fields.

# Configuration File Structure
The TOML configuration file should include sections like:
    - **spins**: Defines spin properties.
    - **grape_parameters**: Parameters for the GRAPE optimization.
    - **optimization_parameters**: Parameters for the hyperparameter optimization.
    - **control_field**: Control field specifications.
    - **save_files**: Settings for saving outputs.
    - **plot**: Plot settings.

# Returns
- Produces and saves results depending on configuration settings, including optimized control fields, cost values, optional Bruker export, and plots.

# Example
```julia
run_grape_optimization("path/to/config.toml")
"""
function run_grape_optimization(config_path::String)
    tm = TOML.parsefile(config_path)
    @info "Configuration:"
    pprintln(tm)

    offsets = collect(-tm["spins"]["offset"]:1:tm["spins"]["offset"])

    # Spin Object
    spins = GrapeMR.Spin(
        tm["spins"]["M0"],
        [s["T1"] for s in tm["spins"]["intrinsics"]],
        [s["T2"] for s in tm["spins"]["intrinsics"]],
        offsets, tm["spins"]["delta_B1"],
        [s["target"] for s in tm["spins"]["intrinsics"]],
        [s["label"] for s in tm["spins"]["intrinsics"]]
    )

    # Grape Parameters
    mask_dict = Dict(k => Bool(v) for (k, v) ∈ tm["grape_parameters"]["fields2optimize"])
    grape_params = GrapeParams(
        tm["grape_parameters"]["time_steps"],
        eval(Symbol(tm["grape_parameters"]["cost_function"])),
        mask_dict
    )

    # Optimization Parameters
    if tm["optimization_parameters"]["hyper_opt"]
        hyper_opt = bohb_hyperopt(spins, grape_params, LinRange(0.01, 0.5, 9), 2187)
        # hyper_opt = random_hyperopt(spins, grape_params, LinRange(0.01, 1.0, 15), range(500, 2000, step = 100))
        Tc, poly_start, poly_degree, max_iter = hyper_opt.minimizer
        opt_params = OptimizationParams(
            poly_start,
            poly_degree,
            Int(ceil(max_iter))
        )
        # Initial RF Pulse Object
        control_field = spline_RF(grape_params.N, Tc, tm["control_field"]["B1ref"])
    else
        opt_params = OptimizationParams(
            tm["optimization_parameters"]["poly_start"],
            tm["optimization_parameters"]["poly_degree"],
            Int(ceil(tm["optimization_parameters"]["max_iter"]))
        )
        # Initial RF Pulse Object
        control_field = spline_RF(grape_params.N, tm["control_field"]["control_time"], tm["control_field"]["B1ref"])
    end

    # Parameters 
    params = Parameters(grape_params, opt_params)

    # Run Optimization
    grape_output = grape(params, control_field, spins)

    # Save data
    if tm["save_files"]["enabled"]
        # Save output data
        if tm["optimization_parameters"]["hyper_opt"]
            experiment_folder = save_grape_data(grape_output; folder_path=tm["save_files"]["folder_path"])
            experiment_folder = save_hyperopt_data(hyper_opt; folder_path=tm["save_files"]["folder_path"])
        else
            experiment_folder = save_grape_data(grape_output; folder_path=tm["save_files"]["folder_path"])
        end
        # Export Bruker data
        if tm["save_files"]["export_bruker"]
            export_bruker(grape_output; folder_path=tm["save_files"]["bruker_folder_path"])
        end
    end

    if tm["plot"]
        # Plots
        display(plot_cost_values(grape_output.cost_values, grape_params))
        display(plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats))
        display(plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control))
        # TODO add if s[t1] > 2 plot(iso[end]) to get time dynamics of the second spin
    end
end



