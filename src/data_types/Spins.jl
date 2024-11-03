abstract type Spins end

struct Spin <: Spins
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    #δ::Vector{Float64}
    B0inho::Float64
    B1inho::Float64
    target::String
    label::String
    Nspins::Float64
end

function Spin(M_ini, T1, T2, B0, B1, targets, labels)
    function vector_spins(args)
        _spins = GrapeMR.Spin[]
        n_spins = length(T1)*length(B0)*length(B1);

        t1, t2, tar, lb = args

        for B0_val ∈ B0
           for B1_val ∈ B1
                _spin = GrapeMR.Spin(M_ini, t1, t2, B0_val, B1_val, tar, lb, n_spins)
                push!(_spins, _spin)
            end
        end
        return _spins
    end

    return mapreduce(vector_spins, vcat, zip(T1, T2, targets, labels)) 
end

struct SpinRange <: Spins
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    B0inho::Float64
    B1inho::Float64
    target::String
    label::String
    Nspins::Float64
end

function SpinRange(M_init, T1, T2, B0inho, B1inho, target, label)
    function vector_spins(args)
        spins = GrapeMR.SpinRange[]
        n_spins = length(T1)*length(T2)

        B1inho, tar, lab = args
        for t2 ∈ T2
            for t1 ∈ T1
                spin = GrapeMR.SpinRange(M_init, t1, t2, B0inho, B1inho, tar, lab, n_spins)
                push!(spins, spin)
            end
        end
        return spins
    end

    spins = vcat(map(vector_spins, zip(B1inho, target, label)))
    return spins
end

struct Magnetization{T<:Real, M<:Union{AbstractVector{T}, AbstractMatrix{T}}}
    dynamics::M
end

struct Isochromat{S<:Spins}
    magnetization::Magnetization
    spin::S
end
