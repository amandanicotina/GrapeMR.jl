"""
    bloch_matrix(cf::InitialControlFields, s::Spins)

bloch_matrix
    # Input  
    - cf: (::InitialControlFields) - Control fields struct
    - s:  (::Spins) - Spin struct

    # Output
    - Calculated 4x4 Bloch matrix
"""
function bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, T1::Float64, T2::Float64)
    # □ Make different calculations for different units
    bloch_matrix = 
        [0.0    0.0     0.0     0.0;
         0.0   -2π/T2  -Bz     -B1y;
         0.0   -Bz     -2π/T2   B1x;
         2π/T1  B1y    -B1x    -2π/T1] 
    
    return bloch_matrix
end



"""
forward_propagation
    # Input  
    - cf: (::InitialControlFields) - Control fields struct
    - s:  (::Spins) - Spin struct

    # Output
    - Magnetization vector 4xN
"""
function forward_propagation(cf::InitialControlFields, s::Spins)
    γ = γ_¹H
    Δt_arr  = range(0.0, cf.t_control, length=cf.N+1)
    M       = zeros(Float64, 4, cf.N+1)
    M[:, 1] = [1.0; s.M_init[1]; s.M_init[2]; s.M_init[3]];

    Bz = 0.0
    Bx = cf.B1x_init_control
    By = cf.B1y_init_control

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz, s.T1, s.T2)
        M[:, i+1] = expv(Δt, b_m, M[:, i])
    end

    return M    
end

"""
backward_propagation
    # Input  
    - cf: (::InitialControlFields) - Control fields struct
    - s:  (::Spins) - Spin struct
    - iso: (::Magnetization) - magnetization vector 4xN

    # Output
    - Adjoint state 4xN
"""

function backward_propagation(cf::InitialControlFields, s::Spins, iso::Magnetization)
    t_arr      = range(cf.t_control, 0.0, length=cf.N+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    χ          = zeros(Float64, 4, cf.N+1)
    χ[:, end]    = cost_gradients["Target One Spin"](iso);

    Bz = 0.0
    Bx = cf.B1x_init_control
    By = cf.B1y_init_control

    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz, s.T1, s.T2)
        bloch_matrix_adjoint = adjoint(b_m)
        χ[:, i] = expv(-Δt[i], bloch_matrix_adjoint, χ[:, i+1]) 
    end

    return round.(χ, digits = 5)
end
