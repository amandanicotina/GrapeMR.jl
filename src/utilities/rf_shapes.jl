"""
    spline_RF(N::Int, t_c::Float64)

Generates a cubic spline pulse 

# Arguments
- `N::Int`: Points
- `t_c::Float64`: Shaped pulse time in seconds
- `B1_ref::Float64`: Reference, normally pulse maximum amplitude

# Output
- ControlField struct 

# Example
```julia
bSSFP(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""
function spline_RF(N, t_c, B1ref)
    len = 10
    time = range(0.0, t_c, length=len)
    # B1x
    field_vals = rand(Float64, len)*B1ref
    spline = CubicSpline(time, field_vals)
    t_values = range(0.0, t_c, length=N)
    spline_vec = [spline(ti) for ti in t_values]
    B1x = spline_vec'
    # B1y
    field_vals = rand(Float64, len)*B1ref
    spline = CubicSpline(time, field_vals)
    t_values = range(0.0, t_c, length=N)
    spline_vec = [spline(ti) for ti in t_values]
    B1y = spline_vec'
    # Bz
    Bz = zeros(1, N)
    
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)

    # Spline plot 
    # plot(time, spline[time])
    # scatter!(time, field_vals)  
    # plot!(t_values, spline_vec)
   
    return control_field
end

# function spline_RF(N::Int64, t_c::Float64, B1ref::Float64; BW_Hz=10.0)
#     len = 10
#     time = range(0.0, t_c, length=len)

#     # Smoothing factor based on the desired bandwidth
#     smoothing_factor = BW_Hz / (1.0 / t_c)  # Approximate relationship between BW and smoothing

#     # Generate smoother random field values by controlling variation
#     field_vals = cumsum(randn(len))  # Generate cumulative sum of random numbers
    
#     # Manually calculate mean and standard deviation
#     mean_val = sum(field_vals) / length(field_vals)
#     std_val = sqrt(sum((field_vals .- mean_val).^2) / length(field_vals))

#     # Normalize
#     field_vals = (field_vals .- mean_val) ./ std_val
#     field_vals = field_vals .* B1ref .* smoothing_factor  # Apply B1ref scaling and smoothing

#     # B1x
#     spline = CubicSpline(time, field_vals)
#     t_values = range(0.0, t_c, length=N)
#     spline_vec_x = [spline(ti) for ti in t_values]
#     B1x = spline_vec_x'

#     # B1y (Generate independently, similar to B1x)
#     field_vals = cumsum(randn(len))  # Same idea: generate smooth random values
    
#     # Manually calculate mean and standard deviation for B1y
#     mean_val = sum(field_vals) / length(field_vals)
#     std_val = sqrt(sum((field_vals .- mean_val).^2) / length(field_vals))

#     # Normalize
#     field_vals = (field_vals .- mean_val) ./ std_val
#     field_vals = field_vals .* B1ref .* smoothing_factor

#     spline = CubicSpline(time, field_vals)
#     spline_vec_y = [spline(ti) for ti in t_values]
#     B1y = spline_vec_y'

#     # Bz remains zero
#     Bz = zeros(1, N)

#     control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)
#     return control_field
# end

"""
    hard_RF(N::Int, t_c::Float64)

Generates a cubic spline pulse 

# Arguments
- `N::Int`: Points
- `t_c::Float64`: Shaped pulse time in seconds

# Output
- 1xN array with spline pulse amplitudes

# Example
```julia
bSSFP(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""
function hard_RF(N, t_c, B1ref)
    # B1x
    B1x = B1ref*ones(1, N)
    # B1y
    B1y = B1ref*zeros(1, N)
    # Bz
    Bz = zeros(1, N)
    
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)
   
    return control_field
end


"""
    sinc_RF(N::Int, t_c::Float64, BW_Hz::Real, flip_angle::Float64)

Generates a sinc pulse with bandwidth BW and flip angle flip_angle in radians

# Arguments
- `N::Int`: Points
- `t_c::Float64`: Shaped pulse time in seconds
- `BW_Hz::Real`: Pulse bandwidth
- `flip_angle::Float64`: Flip angle

# Output
- 1xN array with sinc pulse amplitudes

# Example
```julia
bSSFP(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""
function sinc_RF(N::Int, t_c::Float64, B1ref::Float64; BW_Hz = 500.0, α = π/2)
    t_array = range(0.0, stop=t_c, length=N)
    t = t_array .- t_c / 2
    rot = rad2deg(α) / 360
    flip = rot / diff(t)[1]
    x = BW_Hz .* t
    
    # Generate the B1x component using a sinc function
    B1x = ((flip .* sinc.(x)) ./ 2π)'

    # Generate the B1y component with a 90-degree phase shift (Hilbert transform)
    B1y = ((flip .* sinc.(x .+ π/2)) ./ 2π)'
    Bz = zeros(N)

    # Return the ControlField struct with non-zero B1y
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)
    return control_field
end





"""
    bSSFP_RF(N::Int, nTR::Int, α::Real, TR::Float64)

Generates a bSSFP pulse sequence with bandwidth BW and flip angle flip_angle in radians

# Arguments
- `N::Int`: Points
- `nTR::Int`: How many TRs
- `α::Real`: Flip angle
- `TR::Float64: Repetition time TR in seconds

# Output
- Nx1 array with bSSFP pulse amplitudes

# Example
```julia
bSSFP(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""
function bSSFP_RF(N::Int, nTR::Int, α::Real, TR::Float64)
    Δt  = TR/N
    rf0 = (α/2) / (2π*Δt)
    rf  = α / (2π*Δt)
    bSSFP_vec = Float64[]  

    for n ∈ 1:nTR
        if n == 1
            append!(bSSFP_vec, rf0)
            append!(bSSFP_vec, zeros(N))
        elseif n > 1
            append!(bSSFP_vec, rf)
            append!(bSSFP_vec, zeros(N))
        else
            continue
        end
    end

    return [0.0; bSSFP_vec]
end


