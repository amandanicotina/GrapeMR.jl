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


function save_grape_data(gp::GrapeMR.GrapeOutput; folder_path = pwd())
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
    tar = gp.isochromats[1].spin.target
    lab = gp.isochromats[1].spin.label
    ΔB0 = string(abs(Int(round(gp.isochromats[1].spin.B0inho, digits = 1))))
    opt_folder = "$tar _$lab _$ΔB0"
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
        (df_cost, df_control, df_spins) = save_data_frame(gp)
        CSV.write(joinpath(file_path, "dict_cost_values.csv"), df_cost)
        CSV.write(joinpath(file_path, "dict_control_field.csv"), df_control)
        CSV.write(joinpath(file_path, "dict_iso_spins.csv"), df_spins)
    catch
        error("Unable to save data to CSV files.")
        return
    end

end

function save_data_frame(gp::GrapeMR.GrapeOutput)
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
                                "Chem Shift"   => iso_spin.spin.δ,
                                "Offset [Hz]"  => iso_spin.spin.B0inho,
                                "B1 inho [%]"  => iso_spin.spin.B1inho,
                                "Label"        => iso_spin.spin.label,
                                "Target"       => iso_spin.spin.target,
                                "N° of Spins"  => iso_spin.spin.Nspins) 
                                for iso_spin in gp.isochromats
                        ]      
    # DataFrame format 
    df_cost    = DataFrame(dict_cost_values);
    df_control = DataFrame(dict_control_field)
    df_spins   = DataFrame(dict_iso_spins)
    return df_cost, df_control, df_spins
end
