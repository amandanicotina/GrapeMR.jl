
struct Spin
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    δ::Float64
    B0inho::Float64
    target::String
    label::String
end

struct Magnetization
    # TODO: This could be an NTuple & leverage 4xN known dimension?
    dynamics::Array{Float64}
end

struct Isochromat
    magnetization::Magnetization
    spin::Spin
end
