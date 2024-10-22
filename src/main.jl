using GrapeMR
using TOML

tm = TOML.parsefile("src/default_config.toml")

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
mask = tm["grape_parameters"]["fields2optimize"]
grape_params = GrapeParams(
                tm["grape_parameters"]["time_steps"], 
                Symbol(tm["grape_parameters"]["cost_function"]), 
                reshape(mask, 1, length(mask))  # we need a 1x3 Matrix{Bool} instead of a Vector{Bool}
                )

# Optimization Parameters
if tm["optimization_parameters"]["hyper_opt"]
    hyper_opt = @time bohb_hyperopt(spins, grape_params, LinRange(0.01, 0.5, 9), 2187)
    # random_opt = random_sampler(spins, grape_params, round.(LinRange(0.01, 1.0, 15), digits=2), range(500, 2000, step = 100))
    Tc, poly_start, poly_degree, max_iter = bohb.minimizer 
    opt_params   = OptimizationParams(
                    tm["optimization_parameters"]["poly_start"], 
                    tm["optimization_parameters"]["poly_degree"], 
                    Int(ceil(tm["optimization_parameters"]["max_iter"]))
                    )

else
    opt_params   = OptimizationParams(
                    tm["optimization_parameters"]["poly_start"], 
                    tm["optimization_parameters"]["poly_degree"], 
                    Int(ceil(tm["optimization_parameters"]["max_iter"]))
                    )
end

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse Object
control_field = spline_RF(grape_params.N, tm["control_field"]["control_time"], tm["control_field"]["B1ref"]) 


# Run Optimization
grape_output = @time grape(params, control_field, spins); 

# Save data
if tm["save_files"]["enabled"]
    # Save output data
    if tm["optimization_parameters"]["hyper_opt"]
        experiment_folder = save_output_data(grape_output, hyper_opt; folder_path = tm["save_files"]["folder_path"])
    else
        experiment_folder = save_output_data(grape_output; folder_path = tm["save_files"]["folder_path"])
    end
    # Export Bruker data
    if tm["save_files"]["export_bruker"]
        export_bruker(grape_output; folder_path = tm["save_files"]["bruker_folder_path"])
    end
end

if tm["plot"]
    # Plots
    plot_cost_values(grape_output.cost_values, grape_params)
    plot_magnetization_2D(grape_output.isochromats)
    plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
    plot_control_fields(grape_output.control_field; unit="Hz")
    plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control)
    plot_transverse_time(grape_output.isochromats, grape_output.control_field.t_control)
    plot_transverse_magnetization(grape_output.isochromats)
end


# go = load_grape_data(experiment_folder)
spins = GrapeMR.Spin(
                tm["spins"]["M0"],
                [1e8, 1e8], 
                [1e8, 1e8],
                offsets, tm["spins"]["delta_B1"], 
                [s["target"] for s in tm["spins"]["intrinsics"]],
                [s["label"] for s in tm["spins"]["intrinsics"]]
            )

cf = grape_output.control_field

iso_noRelax = dynamics.(cf, spins)    
plot_magnetization_control_field(grape_output.control_field, iso_noRelax)
plot_transverse_time(iso_noRelax, grape_output.control_field.t_control)
plot_transverse_magnetization(iso_noRelax)


