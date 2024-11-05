"""
    Spins

Abstract type representing a generic spin system in NMR/MRI.
"""
abstract type Spins end

"""
    Spin <: Spins

Gets relaxation parameters and inhomogeneity effects of the spin system.

# Fields
- `M_init::Vector{Float64}`: Initial magnetization vector.
- `T1::Float64`: Longitudinal relaxation time.
- `T2::Float64`: Transverse relaxation time.
- `B0inho::Float64`: B0 inhomogeneity.
- `B1inho::Float64`: B1 inhomogeneity.
- `target::String`: Target magnetization state.
- `label::String`: Label for identifying the spin.
- `Nspins::Float64`: Number of spins in this configuration.
"""
struct Spin <: Spins
    M_init::Vector{Float64}
    T1::Float64
    T2::Float64
    B0inho::Float64
    B1inho::Float64
    target::String
    label::String
    Nspins::Float64
end

"""
    Spin(M_ini, T1, T2, B0, B1, targets, labels)

Constructs a collection of `Spin` instances with specified parameters.

# Arguments
- `M_ini::Vector{Float64}`: Initial magnetization vector.
- `T1::Vector{Float64}`: Longitudinal relaxation times.
- `T2::Vector{Float64}`: Transverse relaxation times.
- `B0::Vector{Float64}`: Array of B0 inhomogeneity values.
- `B1::Vector{Float64}`: Array of B1 inhomogeneity values.
- `targets::Vector{String}`: Target states for each spin configuration.
- `labels::Vector{String}`: Labels for each spin configuration.

# Returns
- A vector of `Spin` instances covering all combinations of the provided parameters.
"""
function Spin(M_ini, T1, T2, B0, B1, targets, labels)
    function vector_spins(args)
        _spins = GrapeMR.Spin[]
        n_spins = length(T1) * length(B0) * length(B1)

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

"""
    Magnetization{T, M}

Represents the magnetization dynamics of a spin system.

# Fields
- `dynamics::M`: Time-evolution data for magnetization, either as a 4xN matrix.
"""
struct Magnetization{T<:Real,M<:Union{AbstractVector{T},AbstractMatrix{T}}}
    dynamics::M
end

"""
    Isochromat{S}

Represents an isochromat, combining magnetization data with a specific spin configuration.

# Fields
- `magnetization::Magnetization`: Magnetization data associated with a collection or single isochromat.
- `spin::S`: Spin configuration for a collection or single isochromat.
"""
struct Isochromat{S<:Spins}
    magnetization::Magnetization
    spin::S
end
