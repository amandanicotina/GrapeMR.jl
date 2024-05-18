function initial_field_spline(N, t_c)
    len = 10
    B1_max = 10;#1/t_c;
    time = range(0.0, t_c, length=len)
    field_vals = rand(Float64, len) .* B1_max
    spline = CubicSpline(time, field_vals)

    t_values = range(0.0, t_c, length=N)
    spline_vec = [spline(ti) for ti in t_values]

    # Spline plot 
    plot(time, spline[time])
    scatter!(time, field_vals)  
    plot!(t_values, spline_vec)
   
    return spline_vec
end



# Make one for sinc



"""
    bSSFP(gp::GrapeMR.GrapeOutput; folder_path = pwd())

Save data related to Grape optimization into files organized in folders.

# Arguments
- `gp::GrapeMR.GrapeOutput`: Grape optimization output.
- `folder_path::String = pwd()`: Folder path where data will be saved.

# Example
```julia
bSSFP(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""

function bSSFP(s::GrapeMR.SteadyState)
    # Spin object (parameters in miliseconds)
    spin = BlochSim.Spin(s.M_init, s.T1*1e3, s.T2*1e3, s.B0inho)

    # RF pulse
    nTR = ceil(3*s.T1/s.TR)
    rf = InstantaneousRF(s.α)
    TR = range(tr+length)

    
end