abstract type AbstractControlField end

"""
    ControlField{T, M1, Mz}

Represents an RF control field with physical units: seconds and Hz.

# Type Parameters
- `T`: Element type (e.g., `Float64`)
- `M1`: Matrix type for the transverse components
- `Mz`: Matrix type for the longitudinal field

# Fields
- `B1x::M1`: x-component of the RF field (in Hz)
- `B1y::M1`: y-component of the RF field (in Hz)
- `B1_ref::T`: Reference RF amplitude (in Hz)
- `Bz::Mz`: Longitudinal field (typically zero, in Hz)
- `t_control::T`: Total pulse duration (in seconds)
"""
mutable struct ControlField{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}} <: AbstractControlField
    B1x::M1
    B1y::M1
    B1_ref::T
    Bz::Mz
    t_control::T
end

"""
    NormalizedControlField{T, M1, Mz}

Represents a unitless RF control field.

Amplitude is normalized by the reference value `B1ref` (in Hz), and time is scaled as
`t ⋅ B1ref`, resulting in a dimensionless representation of the pulse.

Used internally for optimization and gradient calculations.

# Fields
- `B1x::M1`: x-component of the normalized RF field
- `B1y::M1`: y-component of the normalized RF field
- `B1_ref::T`: Always 1.0 (unitless reference)
- `Bz::Mz`: Normalized longitudinal field (typically zero)
- `t_control::T`: Dimensionless total control duration, computed as `t_c * B1ref`
"""
struct NormalizedControlField{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}} <: AbstractControlField
    B1x::M1
    B1y::M1
    B1_ref::T
    Bz::Mz
    t_control::T
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
    normalize_control_field(cf::ControlField) -> NormalizedControlField

Converts a `ControlField` into a unitless `NormalizedControlField`.

- RF amplitude is normalized by `cf.B1_ref` (in Hz)
- Time is normalized as `t ⋅ B1ref`, resulting in a dimensionless control duration

# Arguments
- `cf::ControlField`: Control field in physical units (Hz and seconds)

# Returns
- `NormalizedControlField`: Unitless version of the input field
"""
function normalize_control_field(cf::ControlField)
    B1x_norm = cf.B1x ./ cf.B1_ref
    B1y_norm = cf.B1y ./ cf.B1_ref
    Bz_norm  = cf.Bz  ./ cf.B1_ref
    t_c_norm = cf.t_control * cf.B1_ref
    B1_ref_norm = 1.0  # by definition

    # # DEBUG
    # @info "normalize_control_field" cf.B1_ref maximum_B1x=maximum(abs, cf.B1x) maximum_B1x_norm=maximum(abs, B1x_norm)

    return NormalizedControlField(
        reshape(B1x_norm, 1, :),
        reshape(B1y_norm, 1, :),
        B1_ref_norm,
        reshape(Bz_norm, 1, :),
        t_c_norm,
    )
end

"""
    denormalize_control_field(cf::NormalizedControlField, B1ref::Float64) -> ControlField

Converts a `NormalizedControlField` back into a physical `ControlField`.

- RF amplitude is scaled by `B1ref` (in Hz)
- Time is scaled as `t / B1ref`

# Arguments
- `cf::NormalizedControlField`: Unitless control field
- `B1ref::Float64`: Reference RF amplitude (in Hz)

# Returns
- `ControlField`: Reconstructed field in Hz and seconds
"""
function denormalize_control_field(cf::NormalizedControlField, B1ref::Float64)
    B1x_phys = cf.B1x .* B1ref
    B1y_phys = cf.B1y .* B1ref
    Bz_phys  = cf.Bz  .* B1ref
    t_c_phys = cf.t_control / B1ref

    return ControlField(
        reshape(B1x_phys, 1, :),
        reshape(B1y_phys, 1, :),
        B1ref,
        reshape(Bz_phys, 1, :),
        t_c_phys
    )
end


Base.show(io::IO, cf::NormalizedControlField) =
    print(io, "NormalizedControlField(t_control = $(cf.t_control), B1_ref = $(cf.B1_ref))")

Base.show(io::IO, cf::ControlField) =
    print(io, "ControlField(t_control = $(cf.t_control) s, B1_ref = $(cf.B1_ref) Hz)")