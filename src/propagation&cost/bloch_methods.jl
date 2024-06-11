
"""
    bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

bloch_matrix
    # Input  
    B1x: (::Float64) - B1x step
    B1y: (::Float64) - B1x step
    Bz:  (::Float64) - B1x step
    Γ1:  (::Float64) - B1x step
    Γ2:  (::Float64) - B1x step

    # Output
    - Calculated 4x4 Bloch matrix
"""
function bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, T1::Float64, T2::Float64)

    bloch_matrix = 
        [0.0   0.0   0.0   0.0;
         0.0  -1/T2  Bz   -B1y;
         0.0  -Bz   -1/T2  B1x;
         1/T1  B1y  -B1x  -1/T1] 
    
    return bloch_matrix
end


"""
forward_propagation
    # Input  
    - cf: (::ControlField) - Control fields struct
    - s:  (::Spins) - Spin struct

    # Output
    - Magnetization vector 4xN
"""
function forward_propagation(cf::ControlField, s::Spins)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    M       = zeros(Float64, 4, length(cf.B1x)+1)
    M[:, 1] = [1.0; s.M_init[1]; s.M_init[2]; s.M_init[3]];
    
    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = cf.Bz .+ B0
    Bx = 2π*B1*cf.B1x
    By = 2π*B1*cf.B1y

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], s.T1, s.T2)
        M[:, i+1] = exp(Δt*b_m)*M[:, i]
    end

    return M    
end

"""
backward_propagation
    # Input  
    - cf: (::ControlField) - Control fields struct
    - iso: (::Isochromat) - Magnetization vector 4xN
    - cost_function (::Function) - Cost Function gradient for adjoint state inital state

    # Output
    - Adjoint state 4xN
"""

function backward_propagation(cf::ControlField, iso::Isochromat, cost_grad::Vector{Float64})
    t_arr      = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    χ          = zeros(Float64, 4, length(cf.B1x)+1)
    χ[:, end]  = cost_grad;
    s          = iso.spin

    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = 2π*cf.Bz .+ B0
    Bx = 2π*B1*cf.B1x
    By = 2π*B1*cf.B1y

    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], s.T1, s.T2)
        bloch_matrix_adjoint = adjoint(b_m)
        χ[:, i] = exp(Δt[i]*bloch_matrix_adjoint)*χ[:, i+1]
    end

    return round.(χ, digits = 5)
end


