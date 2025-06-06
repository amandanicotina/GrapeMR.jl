"""
    ControlField{T, M1, Mz}

Represents the RF control field parameters for an NMR/MRI pulse sequence.

# Type Parameters
- `T`: Numeric type (e.g., `Float64`).
- `M1`: Matrix type for the transverse RF components.
- `Mz`: Matrix type for the longitudinal magnetic field component.

# Fields
- `B1x::M1`: x-component of the transverse RF field (1×N matrix).
- `B1y::M1`: y-component of the transverse RF field (1×N matrix).
- `B1_ref::T`: Reference amplitude of the RF field.
- `Bz::Mz`: Longitudinal field component (typically zeros).
- `t_control::T`: Total control duration in seconds.
"""
mutable struct ControlField{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}}
    B1x::M1
    B1y::M1
    B1_ref::T
    Bz::Mz
    t_control::T
end

################################################################################
#                                Helper Functions                              #
################################################################################

"""
    create_spline(spline_time, control_time_vals, B1_vals; rng=Random.GLOBAL_RNG)

Creates a cubic spline interpolation from the specified control points and
evaluates it at given control time values.

# Arguments
- `spline_time::AbstractVector`: Time vector for knot positions.
- `control_time_vals::AbstractVector`: Time vector to evaluate the spline on.
- `B1_vals::AbstractVector`: Amplitude values for interpolation knots.

# Keywords
- `rng`: Random number generator (default: `Random.GLOBAL_RNG`).

# Returns
- Vector of interpolated values at `control_time_vals`.
"""
function create_spline(
    spline_time::AbstractVector,
    control_time_vals::AbstractVector,
    B1_vals::AbstractVector;
    rng = Random.GLOBAL_RNG,
)
    spline = CubicSpline(spline_time, B1_vals)
    return map(t -> spline(t), control_time_vals)
end

"""
    build_control_field(B1x, B1y, B1ref, t_c; Bz=nothing)

Assembles and returns a `ControlField` struct from RF field vectors.

# Arguments
- `B1x::Vector`: Vector of x RF field values.
- `B1y::Vector`: Vector of y RF field values.
- `B1ref::Float64`: Reference amplitude.
- `t_c::Float64`: Total control duration in seconds.

# Keywords
- `Bz`: Optional longitudinal component. Defaults to zero.

# Returns
- A `ControlField` struct with reshaped field arrays.
"""
function build_control_field(
    B1x::Vector,
    B1y::Vector,
    B1ref::Float64,
    t_c::Float64;
    Bz = nothing
)
    B1x_mat = reshape(B1x, 1, :)
    B1y_mat = reshape(B1y, 1, :)
    Bz_mat  = isnothing(Bz) ? zeros(1, length(B1x)) : reshape(Bz, 1, :)
    return ControlField(B1x_mat, B1y_mat, B1ref, Bz_mat, t_c)
end
