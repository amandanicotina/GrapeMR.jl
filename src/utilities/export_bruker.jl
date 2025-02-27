using Printf

"""
    export_bruker(go::GrapeOutput; folder_path = pwd())

Exports the GRAPE optimization results to a Bruker-compatible `.exc` file format for use in TopSpin software.

# Arguments
- `go::GrapeOutput`: The output struct from GRAPE optimization containing control field data.
- `folder_path::String=pwd()`: The directory path where the `.exc` file will be saved. Defaults to the current working directory.

# File Structure
- The generated file includes metadata such as date, time, amplitude, and phase ranges.
- Amplitude and phase data points are formatted for compatibility with Bruker TopSpin.

# Outputs
- Saves a `.exc` file named based on the GRAPE output and current configuration in the specified folder.
"""
function export_bruker(go::GrapeMR.GrapeOutput; folder_path=pwd()) # Add file_path as keyword argument
    # Get control field data
    (amplitudes, phases) = bruker_normalized_amplitudes_and_phases(go.control_field)

    minimum_amplitude = @sprintf("%e", minimum(amplitudes))
    maximum_amplitude = @sprintf("%e", maximum(amplitudes))
    minimum_phases = @sprintf("%e", minimum(phases))
    maximum_phases = @sprintf("%e", maximum(phases))
    integfac = @sprintf("%e", integral_factor(go.control_field))

    # Get current date and time
    current_datetime = now()
    current_date = Dates.format(current_datetime, "dd/mm/yyyy")
    current_time = Dates.format(current_datetime, "HH:MM:SS")

    # Header with general information
    header = string(
        "##TITLE= Parameter file, TopSpin 4.4.0\n",
        "##JCAMPDX= 5.0\n",
        "##DATATYPE= Parameter Values\n",
        "##ORIGIN= TUM\n",
        "##OWNER= <Amanda Nicotina>\n",
        "##DATE= ", current_date, "\n",
        "##TIME= ", current_time, "\n",
        "##MINX= ", minimum_amplitude, "\n",
        "##MAXX= ", maximum_amplitude, "\n",
        "##MINY= ", minimum_phases, "\n",
        "##MAXY= ", maximum_phases, "\n",
        "##\$SHAPE_AMP= 100\n",
        "##\$SHAPE_BWFAC= 0\n",
        "##\$SHAPE_EXMODE= <Excitation>\n",
        "##\$SHAPE_TOTROT= 0.000000e+01\n",
        "##\$SHAPE_TYPE= <Excitation>\n",
        "##\$SHAPE_INTEGFAC= ", integfac, "\n",
        "##\$SHAPE_REPHFAC= 0\n",
        "##\$SHAPE_LENGTH= ", go.control_field.t_control * 1e6, "\n",
        "##\$SHAPE_MODE= 0\n",
        "##NPOINTS= ", length(go.control_field.B1x), "\n",
        "##XYPOINTS=(XY..XY)\n"
    )

    # Amplitude and phases
    amplitudes_and_phases_list = [(@sprintf "%e" amplitude) * ", " * (@sprintf "%e" phase) for (amplitude, phase) in zip(amplitudes, phases)]

    # File content
    content = header * join(amplitudes_and_phases_list, "\n") * "\n ##END"

    # File path and name
    file_name = file_name_string(go) * ".exc"
    file_path = folder_path * file_name
    println(file_path)


    # Open file and write content
    open(file_path, "w") do file
        write(file, content)
    end

end


# TODO export .dat file
# TODO export .CSV file