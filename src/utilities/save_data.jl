
function save_grape_data(gp::GrapeOutput, bool::Bool; folder_path = "/Users/amandanicotina/Documents/ProgressReports/ResultsGrapeMR")
    if bool
        folder_name = Dates.format(today(), "yy-mm-dd")
        full_folder_path = joinpath(folder_path, folder_name)

        # Create/Check folder
        if !isdir(full_folder_path)
            mkdir(full_folder_path)
        end
        # File name
            tar = gp.isochromats[1].spin.target
            lab = gp.isochromats[1].spin.label
            ΔB0 = string(abs(round(gp.isochromats[1].spin.B0inho, digits = 1)))
        # Save data
        file_name = tar * "_" * lab * "_" * ΔB0 ##### => think of a better name
        file_path = joinpath(full_folder_path, file_name)
        open(file_path, "w") do f
            Serialization.serialize(f, gp)
        end
    end
end

function open_grape_data()
    open(file_path, "r") do f
        deserialized_object = deserialize(f, gp)
        println(deserialized_object)
    end
    
end
