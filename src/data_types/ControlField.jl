abstract type AbstractControlField end

"""
    ControlField{T, M1, Mz}

Represents an RF control field with physical units (SI): seconds, Tesla, etc.

# Type Parameters
- `T`: Element type (e.g., `Float64`)
- `M1`: Matrix type for the transverse components
- `Mz`: Matrix type for the longitudinal field

# Fields
- `B1x::M1`: x-component of the RF field
- `B1y::M1`: y-component of the RF field
- `B1_ref::T`: Reference RF amplitude
- `Bz::Mz`: Longitudinal field (typically zero)
- `t_control::T`: Total pulse duration in seconds
"""
mutable struct ControlField{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}} <: AbstractControlField
    B1x::M1
    B1y::M1
    B1_ref::T
    Bz::Mz
    t_control::T
end

"""
    NormalizedControlField

Represents a unitless RF control field. Time and amplitude are scaled to a.u.
Used internally for optimization and gradient calculations.
"""
struct NormalizedControlField{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}} <: AbstractControlField
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
function build_control_field(B1x::Vector, B1y::Vector, B1ref::Float64, t_c::Float64; Bz=nothing)
    B1x_mat = reshape(B1x, 1, :)
    B1y_mat = reshape(B1y, 1, :)
    Bz_mat = isnothing(Bz) ? zeros(1, length(B1x)) : reshape(Bz, 1, :)
    return ControlField(B1x_mat, B1y_mat, B1ref, Bz_mat, t_c)
end



"""
    normalize_control_field(cf::ControlField; t_unit=1e-3, B1_unit=1.0)

Convert a `ControlField` into a `NormalizedControlField` using reference units.

# Arguments
- `t_unit`: Time unit (default: 1e-3 s = 1 ms)
- `B1_unit`: Amplitude unit (default: 1.0)

# Returns
- `NormalizedControlField`
"""
function normalize_control_field(cf::ControlField; t_unit=1e-3, B1_unit=1.0)
    B1x_norm = reshape(cf.B1x ./ B1_unit, 1, :)
    B1y_norm = reshape(cf.B1y ./ B1_unit, 1, :)
    Bz_norm  = reshape(cf.Bz  ./ B1_unit, 1, :)
    t_c_norm = cf.t_control / t_unit
    B1_ref_norm = cf.B1_ref / B1_unit

    # 👇 DEBUG
    @info "normalize_control_field" B1_unit cf.B1_ref maximum_B1x=maximum(abs, cf.B1x) maximum_B1x_norm=maximum(abs, B1x_norm)

    return NormalizedControlField(B1x_norm, B1y_norm, B1_ref_norm, Bz_norm, t_c_norm)
end


"""
    denormalize_control_field(cf::NormalizedControlField; t_unit=1e-3, B1_unit=1.0)

Converts a `NormalizedControlField` back to SI units in a `ControlField`.

# Arguments
- `t_unit`: Time unit (default: 1e-3 s)
- `B1_unit`: Amplitude unit (default: 1.0)

# Returns
- `ControlField`
"""

"""
    denormalize_control_field(cf::NormalizedControlField; t_unit=1e-3, B1_unit=1.0)

Converts a `NormalizedControlField` back to SI units in a `ControlField`.

# Arguments
- `t_unit`: Time unit (default: 1e-3 s)
- `B1_unit`: Amplitude unit (default: 1.0)

# Returns
- `ControlField`
"""
function denormalize_control_field(cf::NormalizedControlField; t_unit=1e-3, B1_unit=1.0)
    B1x_phys = reshape(cf.B1x .* B1_unit, 1, :)
    B1y_phys = reshape(cf.B1y .* B1_unit, 1, :)
    Bz_phys  = reshape(cf.Bz  .* B1_unit, 1, :)
    t_c_phys = cf.t_control * t_unit
    B1_ref_phys = cf.B1_ref * B1_unit
    return ControlField(B1x_phys, B1y_phys, B1_ref_phys, Bz_phys, t_c_phys)
end

# ------------------------------------------------------------------------------
# Optional show methods for clarity
# ------------------------------------------------------------------------------

Base.show(io::IO, cf::NormalizedControlField) =
    print(io, "NormalizedControlField(t_control = $(cf.t_control), B1_ref = $(cf.B1_ref))")

Base.show(io::IO, cf::ControlField) =
    print(io, "ControlField(t_control = $(cf.t_control), B1_ref = $(cf.B1_ref))")