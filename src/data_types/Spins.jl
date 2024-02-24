
struct Spins
    M_init
    T1::Float64
    T2::Float64
    δ::Float64
    target::String
end

# □ 3D array for different δ's
struct Magnetization
    # TODO: This could be an NTuple & leverage 4xN known dimension
    magnetization::Array{Float64}
    spin::Spins
end
