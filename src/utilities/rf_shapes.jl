function spline_RF(N, t_c)
    len = 10
    time = range(0.0, t_c, length=len)
    field_vals = rand(Float64, len)
    spline = CubicSpline(time, field_vals)

    t_values = range(0.0, t_c, length=N)
    spline_vec = [spline(ti) for ti in t_values]

    # Spline plot 
    # plot(time, spline[time])
    # scatter!(time, field_vals)  
    # plot!(t_values, spline_vec)
   
    return spline_vec
end



function sinc_RF(N::Int, t_c::Float64, BW_Hz::Real, α::Float64)
    time = range(0.0, t_c, N);
    t    = time .- t_c/2;
    rot  = rad2deg(α)/360;
    flip = rot/diff(t)[1];
    x    = BW_Hz.*t;
    B1   = ((flip.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
    return B1
end



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

function bSSFP_RF(N::Int, nTR::Int, α::Float64, TR::Float64)
    Δt  = TR/N
    rf0 = (α/2) / (2π*Δt)
    rf  = α / (2π*Δt)
    bSSFP_vec = Float64[]  

    for n ∈ 1:nTR
        if n == 1
            append!(bSSFP_vec, rf0)
            append!(bSSFP_vec, zeros(N))
            @show bSSFP_vec
        elseif n > 1
            append!(bSSFP_vec, rf)
            append!(bSSFP_vec, zeros(N))
        else
            continue
        end
    end

    return [0.0; bSSFP_vec]
end


