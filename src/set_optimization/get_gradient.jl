
function backward_propagation(cf::InitialControlFields, s::Spins, iso::Magnetization)
    γ = γ_¹H

    t_arr      = range(cf.t_control, 0.0, length=cf.N+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    χ          = zeros(Float64, 4, cf.N+1)
    χ[:, end]    = cost_gradients["Grad Euclidean Norm"](iso);

    Bz = 0.0
    Bx = cf.B1x_init_control
    By = cf.B1y_init_control

    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz, s.T1, s.T2)
        bloch_matrix_adjoint = adjoint(b_m)
        χ[:, i] = expv(-Δt[i], bloch_matrix_adjoint, χ[:, i+1]) 
    end

    return χ
end

function control_field_derivative()
end

function gradient_controls(cf::InitialControlFields, s::Spins, iso::Magnetization)
    γ = γ_¹H
    χ = backward_propagation(cf, s, iso)
    Ix = [0 0 0 0; 0 0 0 0; 0 0 0 1; 0 0 -1 0];
    # Iy = [0 0 0 0; 0 0 0 -1; 0 0 0 0; 0 1 0 0];
    t_arr = range(cf.t_control, 0.0, length=cf.N+1)
    Δt    = diff(t_arr)[1]
    M     = magnetization_ODE(cf, s)
    ΔₓJ = zeros(Float64, 1, cf.N)
    for i ∈ 1:cf.N
        ΔₓJ[1,i] = transpose(χ[:,i+1])*Ix*M[:,i]*Δt
    end
    return ΔₓJ
end

function update_control_field(cf::InitialControlFields, s::Spins, iso::Magnetization, ϵ::Float64)
    #∂J = finite_difference_field(gradient_controls, cf, s, iso, ϵ)
    ΔₓJ = gradient_controls(cf, s, iso)
    Bx = cf.B1x_init_control .+ ϵ*ΔₓJ
    #Bx_up = abs.(cf.B1_initial_control[1,2:end]) .+ ϵ*∂J
    #Bx = [cf.B1_initial_control[1,1]; Bx_up]
    return Bx
end
