const γ_¹H = 42.5774688e6 #[Hz/T] 
abstract type ControlField end

struct InitialControlFields
    # □ Change to just one B1 and get the angle and 
    N::Int64
    B1_initial_control::AbstractArray
    B1x_max_amplitude::Float64
    B1y_max_amplitude::Float64
    t_control::Float64
    band_width::Float64
    band_width_step::Float64
end

mutable struct ControlFields
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
end



# □ 3D array for different δ's
struct Magnetization{N}
    magnetization::NTuple{N, Array{Float64}}
    spin::NTuple{N, Spins}
end









# To-Do
### □ Add struct with export Information