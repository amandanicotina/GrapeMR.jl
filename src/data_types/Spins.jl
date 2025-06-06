"""
    Spins

Abstract type representing a generic spin system in NMR/MRI.
"""
abstract type Spins end

"""
    Spin <: Spins

Represents a spin system with relaxation parameters and inhomogeneities.

# Fields
- `m_init::AbstractVector{<:Real}`: Initial magnetization vector.
- `t1::Float64`: Longitudinal relaxation time.
- `t2::Float64`: Transverse relaxation time.
- `b0_inho::Float64`: B0 inhomogeneity.
- `b1_inho::Float64`: B1 inhomogeneity.
- `target::Symbol`: Target magnetization state.
- `label::Symbol`: Label for identifying the spin.
- `n_spins::Int`: Number of spins in this configuration.
"""
struct Spin <: Spins
    m_init::AbstractVector{<:Real}
    T1::Float64
    T2::Float64
    b0_inho::Float64
    b1_inho::Float64
    target::String
    label::String
    n_spins::Int
end

"""
    generate_spins(m_init, t1s, t2s, b0s, b1s, targets, labels)

Generates a flat `Vector{Spin}` containing all combinations of spin parameters.

# Arguments
- `m_init::AbstractVector{<:Real}`: Initial magnetization vector.
- `t1s::Vector{Float64}`: Longitudinal relaxation times.
- `t2s::Vector{Float64}`: Transverse relaxation times.
- `b0s::Vector{Float64}`: B0 inhomogeneity values.
- `b1s::Vector{Float64}`: B1 inhomogeneity values.
- `targets::Vector{Symbol}`: Target states for each spin configuration.
- `labels::Vector{Symbol}`: Labels for each spin configuration.

# Returns
- `Vector{Spin}`: A flat vector of all generated spin configurations.
"""
function generate_spins(m_init, t1s, t2s, b0s, b1s, targets, labels)
    n_spins = length(t1s) * length(b0s) * length(b1s)
    return vcat((
        Spin(m_init, t1, t2, b0, b1, target, label, n_spins)
        for (t1, t2, target, label) in zip(t1s, t2s, targets, labels),
            b0 in b0s, b1 in b1s
    )...)
end





"""
    Magnetization{T, M}

Represents the magnetization dynamics of a spin system.

# Fields
- `dynamics::M`: Time-evolution data for magnetization, either as a vector or matrix.
"""
struct Magnetization{T<:Real, M<:Union{AbstractVector{T}, AbstractMatrix{T}}}
    dynamics::M
end

"""
    Isochromat{S}

Represents an isochromat, combining magnetization data with a specific spin configuration.

# Fields
- `magnetization::Magnetization`: Magnetization data for the isochromat.
- `spin::S`: Spin configuration associated with the magnetization.
"""
struct Isochromat{S<:Spins}
    magnetization::Magnetization
    spin::S
end
