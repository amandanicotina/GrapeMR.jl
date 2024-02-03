
function backward_propagation(cf::InitialControlFields, s::Spins, iso::Magnetization)
    γ = γ_¹H

    t_arr      = range(cf.t_control, 0.0, length=cf.N+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    χ          = zeros(Float64, 4, cf.N+1)
    χ[:, end]    = cost_gradients["Grad Euclidean Norm"](iso);

    Bz = 0.0
    Bx = abs.(cf.B1_initial_control)
    By = angle.(cf.B1_initial_control)

    for i in back_steps:-1:1 
        bloch_matrix = [0.0    0.0        0.0       0.0;
                        0.0    -1/s.T2   -γ*Bz     -γ*By[i];
                        0.0    -γ*Bz     -1/s.T2    γ*Bx[i];
                        1/s.T1  γ*By[i] -γ*Bx[i] -1/s.T1]

        bloch_matrix_adjoint = adjoint(bloch_matrix)
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
    M = magnetization_ODE(cf, s)
    ∂J = transpose(χ[:,end])*Ix*M[:,end] 

    return ∂J
end

function update_control_field(cf::InitialControlFields, s::Spins, iso::Magnetization, ϵ::Float64)
    #∂J = finite_difference_field(gradient_controls, cf, s, iso, ϵ)
    ∂J = gradient_controls(cf, s, iso)
    Bx = abs.(cf.B1_initial_control) .- ϵ*∂J
    #Bx_up = abs.(cf.B1_initial_control[1,2:end]) .+ ϵ*∂J
    #Bx = [cf.B1_initial_control[1,1]; Bx_up]
    return Bx
end
