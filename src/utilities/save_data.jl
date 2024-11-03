"""
    create_folder(path::String)

Creates a folder at the specified path if it doesn't already exist.

# Arguments
- `path::String`: The path where the folder will be created.

# Output
Returns the path to the created folder. If the folder already exists, it simply returns the path.
"""
function create_folder(path::String)
    try
        if !isdir(path)
            mkdir(path)
        end
    catch
        error("Unable to create folder: $path.")
        return
    end
    return path
end

function create_base_folder(folder_path::String)
    folder_name = Dates.format(today(), "yyyy-mm-dd")
    return create_folder(joinpath(folder_path, folder_name))
end


"""
    save_output_data(go::GrapeMR.GrapeOutput; folder_path = pwd())

Save data related to GRAPE optimization into a folder organized by date.

# Arguments
- `go::GrapeMR.GrapeOutput`: The output from a GRAPE optimization process.
- `folder_path::String = pwd()`: The folder path where data will be saved. Defaults to the current working directory.

# Output
Returns the full path to the folder where the data was saved.

# Saved Files
- `grape_output.jld2`: Contains the entire `GrapeOutput` struct in JLD2 format.
- `dict_cost_values.csv`: CSV file containing the cost values from the optimization process.
- `dict_control_field.csv`: CSV file containing control field values (B1x, B1y, Bz, and RF time).
- `dict_iso_spins.csv`: CSV file containing isochromat spin parameters.

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""
function save_grape_data(go::GrapeOutput; folder_path = pwd())
    full_folder_path = create_base_folder(folder_path)

    # Optimization-specific folder hh-mm
    opt_folder = file_name_string(go)
    file_path = create_folder(joinpath(full_folder_path, opt_folder))

    try
        go_dicts = grape_output_to_dict(go)
        save_grape_output(go, file_path)

        # Save the GRAPE data as CSV files
        (df_cost, df_control, df_spins) = gp_dicts_to_data_frame(go_dicts)
        CSV.write(joinpath(file_path, "dict_cost_values.csv"), df_cost)
        CSV.write(joinpath(file_path, "dict_control_field.csv"), df_control)
        CSV.write(joinpath(file_path, "dict_iso_spins.csv"), df_spins)

        return file_path
    catch 
        error("Unable to save GRAPE data")
    end
    return
end

"""
    save_output_data(bohb::GrapeMR.GrapeOutput; folder_path = pwd())

Save data related to BOHB optimization into a folder organized by date.

# Arguments
- `bohb::GrapeMR.GrapeOutput`: The output from a BOHB optimization process.
- `folder_path::String = pwd()`: The folder path where data will be saved. Defaults to the current working directory.

# Output
Returns the full path to the folder where the data was saved.

# Saved Files
- `bohb.jld2`: Contains BOHB result in JLD2 format.
"""
function save_hyperopt_data(hyper_opt; folder_path = pwd())
    full_folder_path = create_base_folder(folder_path)

    try
        save_hyperopt(hyper_opt, full_folder_path)
        return full_folder_path
    catch
        error("Unable to save Hyper-optimization data.")
    end
    return
end

# function save_output_data(go::GrapeOutput, hyper_opt; folder_path = pwd())
#     full_folder_path = create_base_folder(folder_path)

#     # Save GrapeOutput data
#     opt_folder = file_name_string(go)
#     file_path = create_folder(joinpath(full_folder_path, opt_folder))
    
#     try
#         go_dicts = grape_output_to_dict(go)
#         save_grape_output(go, file_path)

#         # Save the GRAPE data as CSV file
#         (df_cost, df_control, df_spins) = gp_dicts_to_data_frame(go_dicts)
#         CSV.write(joinpath(file_path, "dict_cost_values.csv"), df_cost)
#         CSV.write(joinpath(file_path, "dict_control_field.csv"), df_control)
#         CSV.write(joinpath(file_path, "dict_iso_spins.csv"), df_spins)
#         return file_path
#     catch
#         error("Unable to save GRAPE data.")
#         return
#     end

#     # Save the HyperOpt data
#     try
#         save_hyperopt(hyper_opt, file_path)
#         return file_path
#     catch
#         error("Unable to save Hyper-optimization data.")
#     end
# end


function file_name_string(go::GrapeMR.GrapeOutput)
    l = go.isochromats[1].spin.label
    s = go.isochromats[1].spin
    p = go.params.grape_params
    ΔB0 = string(abs(Int(round(go.isochromats[1].spin.B0inho, digits = 1))))
    return string(l, "_", p.cost_function, "_", s.target, "_", "$ΔB0", "Hz")
end

function grape_output_to_dict(gp::GrapeMR.GrapeOutput)
    # Dictionaries with data
    dict_cost_values   = Dict("Cost Values" => gp.cost_values);
    dict_control_field = Dict("B1x [Hz]"    => vec(gp.control_field.B1x),
                              "B1y [Hz]"    => vec(gp.control_field.B1y),  
                              "Bz[Hz]"      => vec(gp.control_field.Bz),
                              "RF Time [s]" => gp.control_field.t_control
                            );
    dict_iso_spins     = [Dict("Initial State" => [iso_spin.spin.M_init],
                                "T1 [s]"       => iso_spin.spin.T1,
                                "T2 [s]"       => iso_spin.spin.T2,
                                "Offset [Hz]"  => iso_spin.spin.B0inho,
                                "B1 inho [%]"  => iso_spin.spin.B1inho,
                                "Label"        => iso_spin.spin.label,
                                "Target"       => iso_spin.spin.target,
                                "N° of Spins"  => iso_spin.spin.Nspins) 
                                for iso_spin in gp.isochromats
                        ]
    return (dict_cost_values, dict_control_field, dict_iso_spins)
end

function gp_dicts_to_data_frame(gp_dicts::Tuple{Dict, Dict, Vector})
    # Dictionaries with data
    (dict_cost_values, dict_control_field, dict_iso_spins) = gp_dicts
    # DataFrame formatS
    df_cost    = DataFrame(dict_cost_values);
    df_control = DataFrame(dict_control_field)
    df_spins   = DataFrame(dict_iso_spins)
    return df_cost, df_control, df_spins
end

function save_gp_dicts(gp_dicts::Tuple{Dict, Dict, Vector}, file_path::String)
    (dict_cost_values, dict_control_field, dict_iso_spins) = gp_dicts
    gp_output_dict = Dict(
        "cost_values"   => dict_cost_values,
        "control_field" => dict_control_field,
        "iso_spins"     => dict_iso_spins
    )
    save(joinpath(file_path, "grape_output_dict.jld2"), "gp_output_dict", gp_output_dict)
end

function save_grape_output(grape_output::GrapeMR.GrapeOutput, file_path::String)
    JLD2.@save joinpath(file_path, "grape_output.jld2") grape_output
end

function save_hyperopt(hyper_opt, file_path::String)
    JLD2.@save joinpath(file_path, "hyper_opt.jld2") hyper_opt
end

function load_grape_data(folder_path::String)
    JLD2.@load joinpath(folder_path, "grape_output.jld2") grape_output
    return grape_output
end

function load_hyperopt_data(folder_path::String)
    JLD2.@load joinpath(folder_path, "hyper_opt.jld2") hyper_opt
    return hyper_opt
end