
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
