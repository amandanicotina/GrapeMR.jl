
function create_spline(spline_time::AbstractArray, control_time_vals::AbstractArray, B1_random_vals::AbstractArray)
    spline = CubicSpline(spline_time, B1_random_vals)
    return map(control_time -> spline(control_time), control_time_vals)
end

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
    spline_time = range(0.0, t_c, length=len)
    control_time = range(0.0, t_c, length=N)
    B1_random_vals = rand(Float64, len)

    # B1x
    B1x = B1ref*create_spline(spline_time, control_time, B1_random_vals)
    B1x_mat = reshape(B1x, 1, :)

    # B1y
    B1y = B1ref*create_spline(spline_time, control_time, B1_random_vals)
    B1y_mat = reshape(B1y, 1, :)

    # Bz
    Bz = zeros(1, N)

    # plot(control_time, B1x)
    # scatter!(spline_time, B1_random_vals)

    return ControlField(B1x_mat, B1y_mat, B1ref, Bz, t_c)
end


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
    B1x = B1ref * ones(1, N)
    # B1y
    B1y = B1ref * zeros(1, N)
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
function sinc_RF(N::Int, t_c::Float64, B1ref::Float64; α=π / 2)
    BW_Hz = 100.0
    t_array = range(0.0, stop=t_c, length=N)
    t = t_array .- t_c / 2
    rot = rad2deg(α) / 360
    flip = rot / diff(t)[1]
    x = BW_Hz .* t

    # B1x 
    B1x = (flip .* sinc.(x)) ./ 2π
    B1x_mat = reshape(B1x, 1, :)

    # B1y
    B1y = (flip .* sinc.(x .+ π / 2)) ./ 2π
    B1y_mat = reshape(B1y, 1, :)

    # B1z
    Bz = zeros(1, N)

    return ControlField(B1x_mat, B1y_mat, B1ref, Bz, t_c)
end

"""
    guassian_RF(N::Int, t_c::Float64, BW_Hz::Real, flip_angle::Float64)

Generates a Gaussian pulse with bandwidth BW and flip angle flip_angle in radians

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
function gaussian_RF(N::Int, t_c::Float64, B1ref::Float64)
    # B1x
    B1x = B1ref * exp.(-0.5 * ((collect(1:N) .- N / 2) / (N / 10)) .^ 2)
    B1x_mat = reshape(B1x, 1, :)

    # B1y
    B1y = B1ref * exp.(-0.5 * ((collect(1:N) .- N / 2) / (N / 10)) .^ 2)
    B1y_mat = reshape(B1y, 1, :)

    # Bz
    Bz = zeros(1, N)

    return ControlField(B1x_mat, B1y_mat, B1ref, Bz, t_c)

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
    Δt = TR / N
    rf0 = (α / 2) / (2π * Δt)
    rf = α / (2π * Δt)
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


