# TODO add examples to docstrings
# TODO refactor bSSFP_RF to export ControlField struct 

"""
    ControlField{T, M1, Mz}

Represents the RF control field parameters for an NMR/MRI sequence.

# Fields
- `B1x::M1`: Matrix for the x-component of the RF field.
- `B1y::M1`: Matrix for the y-component of the RF field.
- `B1_ref::T`: Reference amplitude of the RF field.
- `Bz::Mz`: Matrix for the z-component of the magnetic field.
- `t_control::T`: Total control time for the sequence.
"""
mutable struct ControlField{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}}
    B1x::M1
    B1y::M1
    B1_ref::T
    Bz::Mz
    t_control::T
end


"""
    create_spline(spline_time::AbstractArray, control_time_vals::AbstractArray, B1_random_vals::AbstractArray)

Creates a cubic spline interpolation for the control field based on the specified time and amplitude values.

# Arguments
- `spline_time::AbstractArray`: Array of time points for the spline.
- `control_time_vals::AbstractArray`: Array of control time points for evaluation.
- `B1_random_vals::AbstractArray`: Array of B1 amplitude values for spline generation.

# Returns
- Array of interpolated control field values at each control time.
"""
function create_spline(spline_time::AbstractArray, control_time_vals::AbstractArray, B1_random_vals::AbstractArray)
    spline = CubicSpline(spline_time, B1_random_vals)
    return map(control_time -> spline(control_time), control_time_vals)
end


"""
    spline_RF(N::Int, t_c::Float64, B1ref::Float64)

Generates a cubic spline-based RF pulse.

# Arguments
- `N::Int`: Number of time points.
- `t_c::Float64`: Duration of the pulse in seconds.
- `B1ref::Float64`: Reference amplitude for scaling the pulse.

# Returns
- A `ControlField` struct with `B1x`, `B1y`, and `Bz` components generated using spline interpolation.
"""
function spline_RF(N, t_c, B1ref)
    len = 10
    spline_time = range(0.0, t_c, length=len)
    control_time = range(0.0, t_c, length=N)
    B1_random_vals = rand(Float64, len)

    # B1x
    B1x = B1ref * create_spline(spline_time, control_time, B1_random_vals)
    B1x_mat = reshape(B1x, 1, :)

    # B1y
    B1y = B1ref * create_spline(spline_time, control_time, B1_random_vals)
    B1y_mat = reshape(B1y, 1, :)

    # Bz
    Bz = zeros(1, N)

    return ControlField(B1x_mat, B1y_mat, B1ref, Bz, t_c)
end


"""
    hard_RF(N::Int, t_c::Float64, B1ref::Float64)

Generates a hard RF pulse with constant amplitude in the x-axis and zero amplitude in the y-axis.

# Arguments
- `N::Int`: Number of time points.
- `t_c::Float64`: Duration of the pulse in seconds.
- `B1ref::Float64`: Amplitude of the RF pulse.

# Returns
- A `ControlField` struct with constant `B1x` and zero `B1y` components.
"""
function hard_RF(N, t_c, B1ref)
    # B1x
    B1x = B1ref * ones(1, N)
    # B1y
    B1y = B1ref * zeros(1, N)
    # Bz
    Bz = zeros(1, N) 

    return ControlField(B1x, B1y, B1ref, Bz, t_c)
end


"""
    sinc_RF(N::Int, t_c::Float64, B1ref::Float64; α=π/2)

Generates a sinc RF pulse with a specified flip angle. Bandwidth hardcoded to 100 Hz. B1x and B1y have a π/2 phase difference

# Arguments
- `N::Int`: Number of time points.
- `t_c::Float64`: Duration of the pulse in seconds.
- `B1ref::Float64`: Reference amplitude for scaling the pulse.
- `α::Float64=π/2`: Flip angle in radians.

# Returns
- A `ControlField` struct with `B1x`, `B1y`, and `Bz` components generated as sinc functions.
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
    gaussian_RF(N::Int, t_c::Float64, B1ref::Float64)

Generates a Gaussian-shaped RF pulse.

# Arguments
- `N::Int`: Number of time points.
- `t_c::Float64`: Duration of the pulse in seconds.
- `B1ref::Float64`: Reference amplitude for scaling the pulse.

# Returns
- A `ControlField` struct with Gaussian-distributed `B1x` and `B1y` components.
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

Generates an RF pulse sequence for bSSFP with specified flip angle and repetition time.

# Arguments
- `N::Int`: Number of points per repetition.
- `nTR::Int`: Number of TR periods (repetitions).
- `α::Real`: Flip angle in radians.
- `TR::Float64`: Repetition time in seconds.

# Returns
- A 1D array containing the bSSFP RF pulse sequence.
"""
# function bSSFP_RF(N::Int, nTR::Int, α::Real, TR::Float64)
#     Δt = TR / N
#     rf0 = (α / 2) / (2π * Δt)
#     rf = α / (2π * Δt)
#     bSSFP_vec = Float64[]

#     t_c = nTR * TR

#     for n ∈ 1:nTR
#         if n == 1
#             append!(bSSFP_vec, rf0)
#             append!(bSSFP_vec, zeros(N))
#         elseif n > 1
#             append!(bSSFP_vec, rf)
#             append!(bSSFP_vec, zeros(N))
#         else
#             continue
#         end
#     end
#     b1 = [0.0; bSSFP_vec]
#     B1x = reshape(b1, 1, :)
#     B1y, Bz = zeros(1, N), zeros(1, N)

#     return ControlField(B1x, B1y, 1.0, Bz, t_c)
# end


function bSSFP_RF(N::Int, nTR::Int, α::Real, TR::Float64)
    # Calculate points per TR period
    points_per_TR = N ÷ nTR  # Integer division to ensure even distribution
    Δt = TR / points_per_TR
    
    # Calculate RF amplitudes
    rf0 = (α / 2) / (2π * Δt)
    rf = α / (2π * Δt)
    
    # Initialize vector with zeros
    bSSFP_vec = zeros(N)
    
    # Fill in RF pulses at the start of each TR
    for n in 1:nTR
        idx = (n-1) * points_per_TR + 1  # Index for start of each TR
        if n == 1
            bSSFP_vec[idx] = rf0
        else
            bSSFP_vec[idx] = rf
        end
    end
    
    # Create control field
    t_c = nTR * TR
    B1x = reshape(bSSFP_vec, 1, :)
    B1y = zeros(1, N)
    Bz = zeros(1, N)
    
    return ControlField(B1x, B1y, 1.0, Bz, t_c)
end
