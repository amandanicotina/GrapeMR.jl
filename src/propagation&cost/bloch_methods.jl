function normalization(M_ini, T1, T2, B0, target, label, t_c, B1x, B1y, Bz)
    # □ Use Unitful to normalize based on units of initial RF field
    
    # Omega reference for the normalization
    ω_ref = all(B1x .== 0.0) ? maximum(B1y) : maximum(B1x)

    # Control Field
    τ  = ω_ref*t_c
    uz = Bz./ω_ref
    ux = B1x./ω_ref
    uy = B1y./ω_ref
    ux_max, uy_max = ω_ref, ω_ref

    init_control_field = ControlField(ux, uy, ux_max, uy_max, τ, uz, uz)
    
    # Spins
    function normalized_spin(t1_t2)
        spins = Spin[]  
        n_spins = length(B0);

        t1, t2, tar, lb = t1_t2
        Γ1 = 1/(ω_ref*t1)
        Γ2 = 1/(ω_ref*t2)
        u0 = B0./ω_ref

        for u0_val in u0
            spin = Spin(M_ini, Γ1, Γ2, 0.0, u0_val, tar, lb, n_spins)
            push!(spins, spin)
        end
        return spins
    end
    spins = vcat(map(normalized_spin, zip(T1, T2, target, label))...) 

    return spins, init_control_field
end




"""
    bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

bloch_matrix
    # Input  
    B1x: (::Float64) - Normalized B1x step
    B1y: (::Float64) - Normalized B1x step
    Bz:  (::Float64) - Normalized B1x step
    Γ1:  (::Float64) - Normalized B1x step
    Γ2:  (::Float64) - Normalized B1x step

    # Output
    - Calculated 4x4 Bloch matrix
"""
function bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

    bloch_matrix = 
        [0.0   0.0   0.0   0.0;
         0.0  -Γ2    Bz   -B1y;
         0.0  -Bz   -Γ2    B1x;
         Γ1    B1y  -B1x  -Γ1] 
    
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
function forward_propagation(cf::ControlField, s::Spin)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    M       = zeros(Float64, 4, length(cf.B1x)+1)
    M[:, 1] = [1.0; s.M_init[1]; s.M_init[2]; s.M_init[3]];
    
    ΔB0 = s.B0inho
    Bz  = cf.band_width .+ ΔB0
    Bx  = cf.B1x
    By  = cf.B1y

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = 2π*bloch_matrix(Bx[i], By[i], Bz[i], s.T1, s.T2)
        M[:, i+1] = expv(Δt, b_m, M[:, i])
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

function backward_propagation(cf::ControlField, iso::Isochromat, cost_function::Function)
    t_arr      = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    χ          = zeros(Float64, 4, length(cf.B1x)+1)
    χ[:, end]  = cost_function(iso);
    s          = iso.spin

    ΔB0 = s.B0inho
    Bz  = cf.band_width .+ ΔB0
    Bx  = cf.B1x
    By  = cf.B1y

    for i in back_steps:-1:1 
        b_m = 2π*bloch_matrix(Bx[i], By[i], Bz[i], s.T1, s.T2)
        bloch_matrix_adjoint = adjoint(b_m)
        χ[:, i] = expv(Δt[i], bloch_matrix_adjoint, χ[:, i+1]) 
    end

    return round.(χ, digits = 5)
end


function steady_state(α::Float64, Γ1::Float64, Γ2::Float64, TR::Float64)
    # Calculating signal from geometric derivation
    ϕ = 2π*Δf*TR # Offset angle
    β = 2*atan(tan(α/2) / cos((ϕ)/2)) # β in rads

    Mxy = M0 / (cot(β/2) + (T1/T2)*tan(β/2))
    Mz  = M0*cot(β/2)
    return Mxy, Mz
end