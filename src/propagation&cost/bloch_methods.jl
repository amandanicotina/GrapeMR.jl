"""
    bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

Calculates the Bloch matrix for spin dynamics.

# Arguments
- `B1x::Float64`: x-component of the B1 field.
- `B1y::Float64`: y-component of the B1 field.
- `Bz::Float64`: z-component of the magnetic field.
- `Γ1::Float64`: Longitudinal relaxation rate.
- `Γ2::Float64`: Transverse relaxation rate.

# Returns
- A 4x4 Bloch matrix based on the given field components and relaxation rates.
"""
function bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

    bloch_matrix =
    SA[0.0 0.0 0.0 0.0;
        0.0 -Γ2  -Bz   B1y;
        0.0  Bz  -Γ2  -B1x;
        Γ1  -B1y  B1x -Γ1]

    return bloch_matrix
end

"""
    forward_propagation!(M::AbstractMatrix, cf::ControlField, s::Spins)

In-place version of forward propagation that updates the provided magnetization matrix.

# Arguments  
- `M::AbstractMatrix`: Preallocated 4×N matrix to store the forward-propagated magnetization.
- `cf::ControlField`: Struct containing control field parameters.
- `s::Spins`: Spin struct with initial magnetization and relaxation parameters.

# Returns
- `M`: Updated magnetization matrix (4×N) with forward propagation results.
"""
function forward_propagation!(M::AbstractMatrix, cf::ControlField, s::Spins)
    Δt_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    B1 = s.b1_inho
    B0 = 2π * s.b0_inho

    Bx = 2π * B1 .* cf.B1x
    By = 2π * B1 .* cf.B1y
    Bz = 2π .* cf.Bz .+ B0

    M[:, 1] .= [1.0, s.m_init[1], s.m_init[2], s.m_init[3]]

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        mul!(
            view(M, :, i + 1),
            exp(Δt * b_m),
            view(M, :, i)
        )
    end

    return M
end


"""
    backward_propagation!(χ::AbstractMatrix, cf::ControlField, iso::Isochromat, cost_grad::AbstractVector)

In-place backward propagation for the adjoint state matrix, calculating gradients for the control fields.

# Arguments  
- `χ::AbstractMatrix`: Preallocated matrix (4×N) to store the adjoint state.
- `cf::ControlField`: Struct containing control field parameters.
- `iso::Isochromat`: Isochromat containing the spin configuration.
- `cost_grad::AbstractVector`: Gradient of the cost function for the initial adjoint state.

# Returns
- `χ`: Updated adjoint state matrix (4×N).
"""
function backward_propagation!(
    χ::AbstractMatrix, cf::ControlField, iso::Isochromat, cost_grad::AbstractVector
)
    t_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    Δt = diff(t_arr)
    s = iso.spin

    B1 = s.b1_inho
    B0 = 2π * s.b0_inho
    Bx = 2π * B1 .* cf.B1x
    By = 2π * B1 .* cf.B1y
    Bz = 2π .* cf.Bz .+ B0

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2

    χ[:, end] .= cost_grad

    for i in (length(Δt)):-1:1
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        mul!(
            view(χ, :, i),
            exp(Δt[i] * adjoint(b_m)),
            view(χ, :, i + 1)
        )
    end

    return χ
end

"""
    dynamics(cf::ControlField, spin::Spins)

Computes spin dynamics under a control field, using in-place forward propagation.

# Arguments
- `cf::ControlField`: Control field affecting the spin dynamics.
- `spin::Spins`: Spin object.

# Returns
- `iso::Isochromat`: Resulting isochromat with magnetization trajectory.
"""
function dynamics(cf::ControlField, spin::Spins)
    mag = zeros(Float64, 4, size(cf.B1x, 2) + 1)  # preallocate
    forward_propagation!(mag, cf, spin)
    dyn = GrapeMR.Magnetization(mag)
    return Isochromat(dyn, spin)
end
