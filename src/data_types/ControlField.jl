abstract type ControlField end


struct ControlFields <: ControlField
    B1x::AbstractArray
    B1y::AbstractArray
    B1x_max_amp::Float64
    B1y_max_amp::Float64
    t_control::Float64
    band_width::AbstractArray
    band_width_step::AbstractArray
end