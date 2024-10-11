"""
    save_grape_data(gp::GrapeMR.GrapeOutput; folder_path = pwd())

Save data related to Grape optimization into files organized in folders.

# Arguments
- `gp::GrapeMR.GrapeOutput`: Grape optimization output.
- `folder_path::String = pwd()`: Folder path where data will be saved.

# Example
```julia
save_grape_data(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""

# Add writing a .txt that had all parameters written down?

function save_grape_data(go::GrapeMR.GrapeOutput; folder_path = pwd())
    # General folder yyyy-mm-dd
    folder_name      = Dates.format(today(), "yyyy-mm-dd")
    full_folder_path = joinpath(folder_path, folder_name)   
    # Create/Check folder
    try
        if !isdir(full_folder_path)
            mkdir(full_folder_path)
        end
    catch
        error("Unable to create general folder.")
        return
    end

    # Optimization especific folder hh-mm
    opt_folder = file_name_string(go)
    file_path  = joinpath(full_folder_path, opt_folder)
    # Create/Check folder
    try
        if !isdir(file_path)
            mkdir(file_path)
        end
    catch
        error("Unable to create optimization-specific folder.")
        return
    end

    try
        go_dicts = grape_output_to_dict(go)
        # save_gp_dicts(go_dicts, file_path)
        save_gp_struct(go, file_path)
        (df_cost, df_control, df_spins) = gp_dicts_to_data_frame(go_dicts)
        CSV.write(joinpath(file_path, "dict_cost_values.csv"), df_cost)
        CSV.write(joinpath(file_path, "dict_control_field.csv"), df_control)
        CSV.write(joinpath(file_path, "dict_iso_spins.csv"), df_spins)
        return file_path
    catch
        error("Unable to save data to CSV files.")
        return
    end

end

function grape_output_to_dict(gp::GrapeMR.GrapeOutput)
    # Dictionaries with data
    dict_cost_values   = Dict("Cost Values"   => gp.cost_values);
    dict_control_field = Dict("B1x [Hz]"      => vec(gp.control_field.B1x),
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
        "cost_values" => dict_cost_values,
        "control_field" => dict_control_field,
        "iso_spins" => dict_iso_spins
    )
    save(joinpath(file_path, "grape_output_dict.jld2"), "gp_output_dict", gp_output_dict)
end

function save_gp_struct(grape_output::GrapeMR.GrapeOutput, file_path::String)
    JLD2.@save joinpath(file_path, "grape_output.jld2") grape_output
end

function file_name_string(go::GrapeMR.GrapeOutput)
    l = go.isochromats[1].spin.label
    s = go.isochromats[1].spin
    p = go.params.grape_params
    ΔB0 = string(abs(Int(round(go.isochromats[1].spin.B0inho, digits = 1))))
    return string(l, "_", p.cost_function, "_", s.target, "_", "$ΔB0", "Hz")
end


function load_grape_data(folder_path::String)
    file_name = "grape_output.jld2"
    JLD2.@load joinpath(folder_path, file_name) grape_output
    return grape_output
end