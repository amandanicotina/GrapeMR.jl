
struct Spins
    M_init
    Γ1::Float64
    Γ2::Float64
    δ::Float64
    target::String
end

# □ 3D array for different δ's
struct Magnetization{N}
    magnetization::NTuple{N, Array{Float64}}
    spin::NTuple{N, Spins}
end

