mutable struct ControlField
    B1x::AbstractArray
    B1y::AbstractArray
    B1_ref::Float64
    Bz::AbstractArray
    t_control::Float64
end

# write function ControlField with max Bxmax? : (b_norm = Bxmax) = maximum(B1x) then test to see if now matches