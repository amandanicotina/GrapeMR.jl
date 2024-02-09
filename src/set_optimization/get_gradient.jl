function normalization(cf::InitialControlFields, s::Spins)
    order_of_magnitude = floor(log10(maximum(abs.(cf.B1x_init_control))))
    B_ref = 10^order_of_magnitude;
    
    # Normalizing control fields
    t_norm       = (B_ref*cf.t_control);
    B1x_norm     = (cf.B1x_init_control)./(B_ref);
    B1y_norm     = (cf.B1y_init_control)./(B_ref);
    B1x_norm_max = (2π*cf.B1x_max_amplitude)/(B_ref);
    B1y_norm_max = (2π*cf.B1y_max_amplitude)/(B_ref);

    # Normalizing relaxation values
    T1_norm = 2π/B_ref*s.T1
    T2_norm = 2π/B_ref*s.T2

    field_normalized = InitialControlFields(cf.N, B1x_norm, B1x_norm_max, B1y_norm, B1y_norm_max, t_norm, cf.band_width, cf.band_width_step)
    spin_normalized  = Spins(s.M_init, T1_norm, T2_norm, s.δ, s.target)
    return field_normalized, spin_normalized
end

function gradient_controls(cf::InitialControlFields, s::Spins, iso::Magnetization)
    γ = γ_¹H
    χ = backward_propagation(cf, s, iso)
    Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
    # Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];
    t_arr = range(cf.t_control, 0.0, length=cf.N+1)
    Δt    = diff(t_arr)[1]
    M     = forward_propagation(cf, s)
    ΔxJ = zeros(Float64, 1, cf.N)
    for i ∈ 1:cf.N
        ΔxJ[1,i] = transpose(χ[:,i+1])*Ix*M[:,i]*Δt
    end
    return ΔxJ
end

function update_control_field(cf::InitialControlFields, s::Spins, iso::Magnetization, ϵ::Float64)
    #∂J = finite_difference_field(gradient_controls, cf, s, iso, ϵ)
    ΔxJ = gradient_controls(cf, s, iso)
    Bx = cf.B1x_init_control .+ ΔxJ*ϵ
    #Bx_up = abs.(cf.B1_initial_control[1,2:end]) .+ ϵ*∂J
    #Bx = [cf.B1_initial_control[1,1]; Bx_up]
    return Bx
end
