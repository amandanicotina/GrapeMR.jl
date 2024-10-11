struct SteadyState <: Spins
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    B0inho::Float64
    B1inho::Float64
    target::String
    label::String
    Nspins::Float64
    α::Float64
    Δϕ::Float64
    TR::Float64
    TE::Float64
    M_ss::Vector{Float64}
end

function SteadyState(M_init, T1, T2, B0, B1, target, label, α, Δϕ, TR, TE)
    function vector_spins(args)
        _spins = GrapeMR.SteadyState[]  
        n_spins = length(T1)*length(B0)*length(B1);

        t1, t2, tar, lb = args
        for B0_val in B0
            ss = steady_state_matrix(TR, TE, t1, t2, M_init[3], B0_val, α, Δϕ)
            M_ss = [getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)]
            for B1_val in B1
                _spin = GrapeMR.SteadyState(M_init, t1, t2, B0_val, B1_val, tar, lb, n_spins, α, Δϕ, TR, TE, M_ss)
                push!(_spins, _spin)
            end
        end

        return _spins
    end

    spins = vcat(map(vector_spins, zip(T1, T2, target, label))...) 
    return spins
end


struct SteadyStateData
    offset::Vector{Float64}
    ss_matrix::Vector{ComplexF64}
    ss_manual::Vector{ComplexF64}
    ss_geometric::Vector{Float64}
end