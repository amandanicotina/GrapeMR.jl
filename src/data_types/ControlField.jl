mutable struct ControlField{T<:Real,M1<:AbstractMatrix{T},Mz<:AbstractMatrix{T}}
    B1x::M1
    B1y::M1
    B1_ref::T
    Bz::Mz
    t_control::T
end
