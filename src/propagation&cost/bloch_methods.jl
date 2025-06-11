using StaticArrays

"""
    bloch_matrix(B1x, B1y, Bz, Γ1, Γ2)

Constructs the Bloch matrix for spin dynamics in the rotating frame.

# Arguments
- `B1x`, `B1y`, `Bz`: Magnetic field components (in rad/s).
- `Γ1`: Longitudinal relaxation rate (1/T1).
- `Γ2`: Transverse relaxation rate (1/T2).

# Returns
- A 4×4 matrix that evolves the magnetization vector via `exp(Δt * A)`.
"""
function bloch_matrix(B1x::T, B1y::T, Bz::T, Γ1::T, Γ2::T) where {T<:Real}
    return SA[0.0   0.0   0.0   0.0;
              0.0  -Γ2  -Bz   B1y;
              0.0   Bz  -Γ2  -B1x;
              Γ1  -B1y  B1x  -Γ1]
end


################################################################################
#                             Forward Propagation                              #
################################################################################

"""
    forward_propagation!(M, cf::NormalizedControlField, s)

Simulates forward Bloch dynamics using normalized units.

# Arguments
- `M::Matrix{Float64}`: Preallocated 4×N matrix (magnetization).
- `cf::NormalizedControlField`: Normalized control field.
- `s::Spins`: Spin parameters (T1, T2, initial magnetization, B₀/B₁ inhomogeneity).

# Returns
- The updated magnetization matrix `M`.
"""
function forward_propagation!(M::AbstractMatrix, cf::NormalizedControlField, s::Spins)
    Δt_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    M[:, 1] .= [1.0, s.m_init[1], s.m_init[2], s.m_init[3]]

    Bx = (2π * s.b1_inho) .* cf.B1x
    By = (2π * s.b1_inho) .* cf.B1y
    Bz = 2π .* (cf.Bz .+ s.b0_inho)

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    for (i, Δt) in enumerate(diff(Δt_arr))
        mul!(view(M, :, i + 1), exp(Δt * bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)), view(M, :, i))
    end
    return M
end

"""
    forward_propagation!(M, cf::ControlField, s)

Simulates forward Bloch dynamics using SI units (Tesla and seconds).

# Arguments
- `M::Matrix{Float64}`: Preallocated 4×N matrix (magnetization).
- `cf::ControlField`: Control field with SI units.
- `s::Spins`: Spin parameters (T1, T2, initial magnetization, B₀/B₁ inhomogeneity).

# Returns
- The updated magnetization matrix `M`.
"""
function forward_propagation!(M::AbstractMatrix, cf::ControlField, s::Spins)
    Δt_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    M[:, 1] .= [1.0, s.m_init[1], s.m_init[2], s.m_init[3]]

    Bx = (2π * s.b1_inho) .* cf.B1x
    By = (2π * s.b1_inho) .* cf.B1y
    Bz = 2π .* (cf.Bz .+ s.b0_inho)

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    for (i, Δt) in enumerate(diff(Δt_arr))
        mul!(view(M, :, i + 1), exp(Δt * bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)), view(M, :, i))
    end
    return M
end


################################################################################
#                          Backward Propagation                                # 
################################################################################

"""
    backward_propagation!(χ, cf::NormalizedControlField, iso, cost_grad)

Simulates backward adjoint dynamics using normalized control fields.

# Arguments
- `χ::Matrix{Float64}`: Preallocated 4×N matrix (adjoint state).
- `cf::NormalizedControlField`: Normalized control field.
- `iso::Isochromat`: Contains spin and forward magnetization.
- `cost_grad::Vector{Float64}`: Gradient of the cost at final time.

# Returns
- The updated adjoint state matrix `χ`.
"""
function backward_propagation!(χ::AbstractMatrix, cf::NormalizedControlField, iso::Isochromat, cost_grad::AbstractVector)
    Δt = diff(range(0.0, cf.t_control, length(cf.B1x) + 1))
    s  = iso.spin
    χ[:, end] .= cost_grad

    Bx = (2π * s.b1_inho) .* cf.B1x
    By = (2π * s.b1_inho) .* cf.B1y
    Bz = 2π .* (cf.Bz .+ s.b0_inho)

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    for i in length(Δt):-1:1
        mul!(view(χ, :, i), exp(Δt[i] * adjoint(bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2))), view(χ, :, i + 1))
    end
    return χ
end

"""
    backward_propagation!(χ, cf::ControlField, iso, cost_grad)

Simulates backward adjoint dynamics using SI-unit control fields.

# Arguments
- `χ::Matrix{Float64}`: Preallocated 4×N matrix (adjoint state).
- `cf::ControlField`: Control field in Tesla/seconds.
- `iso::Isochromat`: Contains spin and forward magnetization.
- `cost_grad::Vector{Float64}`: Gradient of the cost at final time.

# Returns
- The updated adjoint state matrix `χ`.
"""
function backward_propagation!(χ::AbstractMatrix, cf::ControlField, iso::Isochromat, cost_grad::AbstractVector)
    Δt = diff(range(0.0, cf.t_control, length(cf.B1x) + 1))
    s  = iso.spin
    χ[:, end] .= cost_grad

    Bx = (2π * s.b1_inho) .* cf.B1x
    By = (2π * s.b1_inho) .* cf.B1y
    Bz = 2π .* (cf.Bz .+ s.b0_inho)

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    for i in length(Δt):-1:1
        mul!(
            view(χ, :, i), 
            exp(Δt[i] * adjoint(bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2))), 
            view(χ, :, i + 1)
        )
    end
    return χ
end


################################################################################
#                                Dynamics                                      #
################################################################################

"""
    dynamics(cf, spin)

Simulates Bloch dynamics using either SI or normalized control fields.

# Arguments
- `cf::AbstractControlField`: A `ControlField` or `NormalizedControlField`.
- `spin::Spins`: Spin configuration.

# Returns
- An `Isochromat` struct containing the forward magnetization.
"""
function dynamics(cf::AbstractControlField, spin::Spins)
    M = zeros(Float64, 4, size(cf.B1x, 2) + 1)
    forward_propagation!(M, cf, spin)
    return Isochromat(GrapeMR.Magnetization(M), spin)
end
