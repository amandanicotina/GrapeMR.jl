function initial_field_spline(N, t_c)
    len = 10
    B1_max = 10;#1/t_c;
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
    B_ref = all(B1x .== 0.0) ? maximum(B1y) : maximum(B1x)

    # Control Field
    τ  = B_ref*t_c
    uz = Bz./B_ref
    ux = B1x./B_ref
    uy = B1y./B_ref
    norm_control_field = ControlField(ux, uy, B_ref, uz, τ) 

    # Spins
    function normalized_spins(t1_t2)
        spins = Spin[]  
        n_spins = length(B0)*length(B1);

        t1, t2, tar, lb = t1_t2
        Γ1 = 1/(B_ref*t1)
        Γ2 = 1/(B_ref*t2)
        u0 = B0./B_ref
        u1 = B1

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

function rescale(gp::GrapeMR.GrapeOutput)
    grape_output_rescale = GrapeOutput([], deepcopy(gp.control_field), gp.cost_values)
    B_ref = gp.control_field.B1_ref

    # Spins 
    for iso ∈ gp.isochromats
        spin = iso.spin
        mag  = Magnetization(iso.magnetization.dynamics)
        T1 = round(1/(B_ref*spin.T1), digits=3)
        T2 = round(1/(B_ref*spin.T2), digits=3)
        B0 = B_ref*(spin.B0inho)
        B1 = B_ref*(spin.B1inho)
        spin_norm = Spin(spin.M_init, T1, T2, 0.0, B0, B1, spin.target, spin.label, spin.Nspins)
        iso_norm = Isochromat(mag, spin_norm)
        push!(grape_output_rescale.isochromats, iso_norm)
    end

    # Control Fields
    all(gp.control_field.B1x .== 0.0) ? norm_B1x = zeros(1, length(gp.control_field.B1x)) : norm_B1x = B_ref*(gp.control_field.B1x./maximum(gp.control_field.B1x))
    all(gp.control_field.B1y .== 0.0) ? norm_B1y = zeros(1, length(gp.control_field.B1y)) : norm_B1y = B_ref*(gp.control_field.B1y./maximum(gp.control_field.B1y))


    grape_output_rescale.control_field.B1x = norm_B1x
    grape_output_rescale.control_field.B1y = norm_B1y
    grape_output_rescale.control_field.t_control = gp.control_field.t_control/B_ref
    return grape_output_rescale

end