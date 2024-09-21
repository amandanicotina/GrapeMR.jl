mutable struct ControlField
    B1x::AbstractArray
    B1y::AbstractArray
    B1_ref::Float64
    Bz::AbstractArray
    t_control::Float64
end

mutable struct ControlFieldNormalized
    B1x::AbstractArray
    B1y::AbstractArray
    B1_ref::Float64
    Bz::AbstractArray
    t_control::Float64
end