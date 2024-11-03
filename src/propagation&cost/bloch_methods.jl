
"""
    bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

# Arguments
- B1x: (::Float64) - B1x step
- B1y: (::Float64) - B1x step
- Bz:  (::Float64) - B1x step
- Γ1:  (::Float64) - B1x step
- Γ2:  (::Float64) - B1x step

# Outputs
- Calculated 4x4 Bloch matrix
"""
function bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

    bloch_matrix = 
        SA[0.0   0.0   0.0   0.0;
           0.0  -Γ2    Bz   -B1y;
           0.0  -Bz   -Γ2   B1x;
           Γ1    B1y  -B1x  -Γ1] 
    
    return bloch_matrix
end

function calculate_bloch_matrix()
    
end

"""
forward_propagation

# Arguments  
- cf: (::ControlField) - Control fields struct
- s:  (::Spins) - Spin struct

# Outputs
- Magnetization vector 4xN
"""
function forward_propagation(cf::ControlField, s::Spins)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    M       = zeros(Float64, 4, length(cf.B1x)+1)
    M[:, 1] = [1.0, s.M_init[1], s.M_init[2], s.M_init[3]]
    
    B0 = 2π.*s.B0inho
    B1 = s.B1inho
    Bz = 2π.*cf.Bz .+ B0
    Bx = (2π*B1).*cf.B1x
    By = (2π*B1).*cf.B1y

    Γ1 = 1/s.T1
    Γ2 = 1/s.T2
    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        mul!(
            view(M, :, i+1), 
            exp(Δt*b_m),
            view(M, :, i)
        )
    end

    return M    
end

function forward_propagation!(M::AbstractMatrix, cf::ControlField, s::Spins)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    
    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = cf.Bz .+ B0
    Bx = 2π*B1*cf.B1x
    By = 2π*B1*cf.B1y

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], s.T1,  s.T2)
        mul!(
            view(M, :, i+1), 
            exp(Δt*b_m),
            view(M, :, i)
        )
    end

    return M    
end

"""
backward_propagation

# Arguments  
- cf: (::ControlField) - Control fields struct
- iso: (::Isochromat) - Magnetization vector 4xN
- cost_function (::Function) - Cost Function gradient for adjoint state inital state

# Outputs
- Adjoint state 4xN
"""

function backward_propagation(cost_grad::AbstractVector, cf::ControlField, iso::Isochromat)
    t_arr      = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    # TODO: refactor this as backward_propagation!(χ, cf, iso, cost_grad)
    χ          = zeros(Float64, 4, length(cf.B1x)+1)
    χ[:, end]  = cost_grad;
    s          = iso.spin
    
    B0 = 2π.*s.B0inho
    B1 = s.B1inho
    Bz = 2π.*cf.Bz .+ B0
    Bx = (2π*B1).*cf.B1x
    By = (2π*B1).*cf.B1y
    
    Γ1 = 1/s.T1
    Γ2 = 1/s.T2
    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        bloch_matrix_adjoint = adjoint(b_m)
        mul!(
            view(χ, :, i),
            exp(Δt[i]*bloch_matrix_adjoint),
            view(χ, :, i+1)
        )
    end

    return round.(χ, digits = 5)
end

function backward_propagation!(χ::AbstractMatrix, cf::ControlField, iso::Isochromat)
    t_arr      = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    s          = iso.spin

    B0 .= 2π.*s.B0inho
    B1 .= s.B1inho
    Bz .= 2π.*cf.Bz .+ B0
    Bx .= (2π*B1).*cf.B1x
    By .= (2π*B1).*cf.B1y
    
    Γ1 = 1/s.T1
    Γ2 = 1/s.T2
    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], Γ1, Γ2)
        bloch_matrix_adjoint = adjoint(b_m)
        mul!(
            view(χ, :, i),
            exp(Δt[i]*bloch_matrix_adjoint),
            view(χ, :, i+1)
        )
    end

    return round.(χ, digits = 5)
end
