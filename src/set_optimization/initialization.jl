function initial_field_spline(N, t_c, B1_max)
    len = 10
    time = range(0.0, t_c, length=len)
    field_vals = rand(Float64, len) .* B1_max
    spline = CubicSpline(time, field_vals)

    t_values = range(0.0, t_c, length=N)
    spline_vec = [spline(ti) for ti in t_values]

    # Spline plot 
    plot(time, spline[time])
    scatter!(time, field_vals)  
    plot!(t_values, spline_vec)
   
    return spline_vec
end

function normalization(M_ini, T1, T2, B0, target, label, t_c, B1x, B1y, B1, Bz)
    # Omega reference for the normalization
    ω_ref = all(B1x .== 0.0) ? maximum(B1y) : maximum(B1x)

    # Control Field
    τ  = ω_ref*t_c
    uz = Bz./ω_ref
    ux = B1x./ω_ref
    uy = B1y./ω_ref
    ux_max, uy_max = ω_ref, ω_ref
    norm_control_field = ControlField(ux, uy, ux_max, uy_max, uz, τ) 

    # Spins
    function normalized_spins(t1_t2)
        spins = Spin[]  
        n_spins = length(B0)*length(B1);

        t1, t2, tar, lb = t1_t2
        Γ1 = 1/(ω_ref*t1)
        Γ2 = 1/(ω_ref*t2)
        u0 = B0./ω_ref
        u1 = B1./ω_ref

        for u0_val ∈ u0
            for u1_val ∈ u1
                spin = Spin(M_ini, Γ1, Γ2, 0.0, u0_val, u1_val, tar, lb, n_spins)
                push!(spins, spin)
            end
        end
        return spins
    end
    spins = vcat(map(normalized_spins, zip(T1, T2, target, label))...) 

    return spins, norm_control_field
end
