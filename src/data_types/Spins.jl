abstract type Spins end

struct Spin <: Spins
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    δ::Float64
    B0inho::Float64
    B1inho::Float64
    target::String
    label::String
    Nspins::Float64
end

function Spin(M_ini, T1, T2, B0, B1, target, label)
    function vector_spins(args)
        spins = GrapeMR.Spin[]  
        n_spins = length(T1)*length(B0)*length(B1);

        t1, t2, tar, lb = args

        for B0_val ∈ B0
           for B1_val ∈ B1
                spin = GrapeMR.Spin(M_ini, t1, t2, 0.0, B0_val, B1_val, tar, lb, n_spins)
                push!(spins, spin)
            end
        end
        return spins
    end

    spins = vcat(map(vector_spins, zip(T1, T2, target, label))...) 
    return spins
end

struct SteadyState <: Spins
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    δ::Float64
    B0inho::Float64
    B1inho::Float64
    target::String
    label::String
    Nspins::Float64
    α::Float64
    Δϕ::Float64
    TR::Float64
    TE::Float64
end

function SteadyState(M_ini, T1, T2, B0, B1, target, label, α, Δϕ, TR, TE)
    function vector_spins(args)
        spins = GrapeMR.SteadyState[]  
        n_spins = length(T1)*length(B0)*length(B1);

        t1, t2, tar, lb = args

        for B0_val ∈ B0
           for B1_val ∈ B1
                spin = GrapeMR.SteadyState(M_ini, t1, t2, 0.0, B0_val, B1_val, tar, lb, n_spins, α, Δϕ, TR, TE)
                push!(spins, spin)
            end
        end
        return spins
    end

    spins = vcat(map(vector_spins, zip(T1, T2, target, label))...) 
    return spins
end


struct Magnetization
    # TODO: This could be an NTuple & leverage 4xN known dimension?
    dynamics::Array{Float64}
end

struct Isochromat
    magnetization::Magnetization
    spin::Spins
end
