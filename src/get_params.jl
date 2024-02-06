const γ_¹H = 42.5774688e6 #[Hz/T] 
abstract type ControlField end

struct InitialControlFields <: ControlField
    N::Int64
    B1x_init_control::AbstractArray
    B1x_max_amplitude::Float64
    B1y_init_control::AbstractArray
    B1y_max_amplitude::Float64
    t_control::Float64
    band_width::Float64
    band_width_step::Float64
end

struct ControlMatrices <: ControlField
    Bx_mat::Array{Float64, 2}
    By_mat::Array{Float64, 2}
    B0_mat::Array{Float64, 2}
end

mutable struct ControlFields <: ControlField
    B1x::AbstractArray
    B1y::AbstractArray
    t_control::Float64
    cost_vals::Array{Float64, 2}
    band_width::Float64
    band_width_step::Float64
end


struct Spins
    M_init
    T1::Float64
    T2::Float64
    δ::Float64
    target::String
end

struct OptimizationParams
    cost_function::String # □ Make a dict with the avaiable cost functions
    max_iter::Int64
    fields_opt::Array{Bool, 3}
    ϵ::Float64
end



# □ 3D array for different δ's
struct Magnetization{N}
    magnetization::NTuple{N, Array{Float64}}
    spin::NTuple{N, Spins}
end









# To-Do
### □ Add struct with export Information