
function magnetization_excitation!(cf::InitialControlFields, s::Spins)
    B₁ = √(cf.B1x_initial_control^2 + cf.(B1y_initial_control)^2)
    Δt = cf.t_control
    α = γ*B₁*Δt # α = γB₁Δt; γ[MHz/T], B₁[T] and Δt[s]
    M
    
end

function magnetization_precession!(cf::InitialControlFields, s::Spins)
    flip = cf.B1x_initial_control
    
end


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
    γ = γ_¹H 

    bloch_matrix = 
        [0.0   0.0     0.0     0.0;
         0.0  -1/T2   -γ*Bz   -γ*B1y;
         0.0  -γ*Bz   -1/T2    γ*B1x;
         1/T1  γ*B1y  -γ*B1x  -1/T1] 
    
    return bloch_matrix
end



"""
magnetization_ODE
    # Input  
    - cf: (::InitialControlFields) - Control fields struct
    - s:  (::Spins) - Spin struct

    # Output
    - Magnetization vector 4xN
"""
function magnetization_ODE(cf::InitialControlFields, s::Spins)
    γ = γ_¹H
    Δt_arr    = range(0.0, cf.t_control, length=cf.N+1)
    mag_exact = zeros(Float64, 4, cf.N+1)
    mag_exact[:, 1] = [1.0; s.M_init[1]; s.M_init[2]; s.M_init[3]];

    Bz = 0.0
    Bx = cf.B1x_init_control
    By = cf.B1y_init_control

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
    bloch_matrix = [0.0    0.0        0.0       0.0;
                    0.0    -1/s.T2   -γ*Bz     -γ*By[i];
                    0.0    -γ*Bz     -1/s.T2    γ*Bx[i];
                    1/s.T1  γ*By[i]  -γ*Bx[i]  -1/s.T1 ] 
    mag_exact[:, i+1] = expv(Δt, bloch_matrix, mag_exact[:, i])
    end

    return mag_exact    
end