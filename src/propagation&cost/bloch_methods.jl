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
            0.0 -Γ2 Bz -B1y;
            0.0 -Bz -Γ2 B1x;
            Γ1 B1y -B1x -Γ1]

    return bloch_matrix
end

"""
    forward_propagation(cf::ControlField, s::Spins)

Performs forward propagation for the magnetization vector under the influence of control fields.

# Arguments  
- `cf::ControlField`: Struct containing control field parameters.
- `s::Spins`: Spin struct with initial magnetization and relaxation parameters.

# Returns
- `M::Matrix{Float64}`: Magnetization matrix (4xN) at each time step.
"""
function forward_propagation(cf::ControlField, s::Spins)
    Δt_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    M = zeros(Float64, 4, length(cf.B1x) + 1)
    M[:, 1] = [1.0, s.M_init[1], s.M_init[2], s.M_init[3]]

    B0 = 2π .* s.B0inho
    B1 = s.B1inho
    Bz = 2π .* cf.Bz .+ B0
    Bx = (2π * B1) .* cf.B1x
    By = (2π * B1) .* cf.B1y

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
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
    forward_propagation!(M::AbstractMatrix, cf::ControlField, s::Spins)

In-place version of forward propagation that updates the provided magnetization matrix.

# Arguments  
- `M::AbstractMatrix`: Matrix to store the forward-propagated magnetization (4xN).
- `cf::ControlField`: Struct containing control field parameters.
- `s::Spins`: Spin struct with initial magnetization and relaxation parameters.

# Returns
- `M`: Updated magnetization matrix (4xN) with forward propagation results.
"""
function forward_propagation!(M::AbstractMatrix, cf::ControlField, s::Spins)
    # Still not working
    Δt_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)

    B0 = 2π * s.B0inho
    B1 = s.B1inho
    Bz = cf.Bz .+ B0
    Bx = 2π * B1 * cf.B1x
    By = 2π * B1 * cf.B1y

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], s.T1, s.T2)
        mul!(
            view(M, :, i + 1),
            exp(Δt * b_m),
            view(M, :, i)
        )
    end

    return M
end

"""
    backward_propagation(cost_grad::AbstractVector, cf::ControlField, iso::Isochromat)

Performs backward propagation for the adjoint state matrix, calculating gradients for the control fields.

# Arguments  
- `cost_grad::AbstractVector`: Gradient of the cost function for the initial adjoint state.
- `cf::ControlField`: Struct containing control field parameters.
- `iso::Isochromat`: Isochromat containing the magnetization data for backward propagation.

# Returns
- `χ::Matrix{Float64}`: Adjoint state matrix (4xN) after backward propagation.
"""
function backward_propagation(cost_grad::AbstractVector, cf::ControlField, iso::Isochromat)
    t_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    Δt = diff(t_arr)
    back_steps = length(Δt)
    # TODO: refactor this as backward_propagation!(χ, cf, iso, cost_grad)
    χ = zeros(Float64, 4, length(cf.B1x) + 1)
    χ[:, end] = cost_grad
    s = iso.spin

    B0 = 2π .* s.B0inho
    B1 = s.B1inho
    Bz = 2π .* cf.Bz .+ B0
    Bx = (2π * B1) .* cf.B1x
    By = (2π * B1) .* cf.B1y

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    for i in back_steps:-1:1
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        bloch_matrix_adjoint = adjoint(b_m)
        mul!(
            view(χ, :, i),
            exp(Δt[i] * bloch_matrix_adjoint),
            view(χ, :, i + 1)
        )
    end

    return round.(χ, digits=5)
end


"""
    backward_propagation!(χ::AbstractMatrix, cf::ControlField, iso::Isochromat)

In-place backward propagation for the adjoint state matrix, calculating gradients for the control fields.

# Arguments  
- `χ::AbstractMatrix`: Matrix to store the adjoint state (4xN).
- `cf::ControlField`: Struct containing control field parameters.
- `iso::Isochromat`: Isochromat containing the magnetization data for backward propagation.

# Returns
- `χ`: Updated adjoint state matrix (4xN) after backward propagation.
"""
function backward_propagation!(χ::AbstractMatrix, cf::ControlField, iso::Isochromat)
    # Still not working
    t_arr = range(0.0, cf.t_control, length(cf.B1x) + 1)
    Δt = diff(t_arr)
    back_steps = length(Δt)
    s = iso.spin

    B0 .= 2π .* s.B0inho
    B1 .= s.B1inho
    Bz .= 2π .* cf.Bz .+ B0
    Bx .= (2π * B1) .* cf.B1x
    By .= (2π * B1) .* cf.B1y

    Γ1 = 1 / s.T1
    Γ2 = 1 / s.T2
    for i in back_steps:-1:1
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        bloch_matrix_adjoint = adjoint(b_m)
        mul!(
            view(χ, :, i),
            exp(Δt[i] * bloch_matrix_adjoint),
            view(χ, :, i + 1)
        )
    end

    return round.(χ, digits=5)
end
