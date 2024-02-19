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
function forward_propagation(cf::ControlFields, s::Spins)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    M       = zeros(Float64, 4, length(cf.B1x)+1)
    M[:, 1] = [1.0; s.M_init[1]; s.M_init[2]; s.M_init[3]];

    Bz = 0.0
    Bx = cf.B1x
    By = cf.B1y

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

function backward_propagation(cf::ControlFields, iso::Magnetization, cost_function::String)
    t_arr      = range(cf.t_control, 0.0, length(cf.B1x)+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    χ          = zeros(Float64, 4, length(cf.B1x)+1)
    χ[:, end]  = cost_gradients[cost_function](iso);
    s = iso.spin[1]

    Bz = 0.0
    Bx = cf.B1x
    By = cf.B1y

    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz, s.T1, s.T2)
        bloch_matrix_adjoint = adjoint(b_m)
        χ[:, i] = expv(-Δt[i], bloch_matrix_adjoint, χ[:, i+1]) 
    end

    return round.(χ, digits = 5)
end
