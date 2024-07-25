"""
Cost Function's Gradients

"""
function cost_function_gradient(iso::Isochromat, cf::Symbol)
    @match cf begin
        :euclidean_norm         => grad_euclidean_norm(iso)
        :target_one_spin        => grad_target_one_spin(iso)
        :target_steady_state    => grad_target_steady_state(iso)
        :targetB0_steady_state  => grad_targetB0_steady_state(iso)
        :target_phase_encoding  => grad_target_phase_encoding(iso)
        :saturation_contrast    => grad_saturation_contrast(iso)
        :saturation_contrast_Mx => grad_saturation_contrast_Mx(iso)
        :saturation_contrast_square       => grad_saturation_contrast_square(iso)
        :saturation_contrast_steady_state => grad_saturation_contrast_steady_state(iso)
    end
    
end
                
function grad_euclidean_norm(iso::Isochromat)
    m = iso.magnetization.dynamics
    s = iso.spin
    norm =  sqrt(sum(m[2:3,end].*m[2:3,end]) + 1e-15)
    Mx_tf = m[2,end]/norm
    My_tf = m[3,end]/norm
    Mz_tf = m[4,end]/norm


    P = [0.0, Mx_tf, My_tf, 0.0]

    return P
end


function grad_target_one_spin(iso::Isochromat; M_tar = [0.0, 1.0, 0.0])  
    m = iso.magnetization.dynamics
    Px  = m[2,end] -  M_tar[1,1]
    Py  = m[3,end] -  M_tar[2,1]
    Pz  = m[4,end] -  M_tar[3,1]

    P_tar = [0.0, Px, Py, Pz]

    return P_tar
end

function grad_target_phase_encoding(iso::Isochromat; M_tar = [0.0, 1.0, 0.0])
    s  = iso.spin
    m  = iso.magnetization.dynamics
    Δt = 0.1 # seconds
    ϕ  = 2π*s.B0inho*Δt

    Px = -sin(ϕ)*(m[2,end] -  M_tar[1,1])/s.Nspins
    Py = -cos(ϕ)*(m[3,end] -  M_tar[2,1])/s.Nspins
    Pz = -(m[4,end] -  M_tar[3,1])/s.Nspins

    P = [0.0, Px, Py, Pz]
    
    return P
end


function grad_saturation_contrast(iso::Isochromat)
    m = iso.magnetization.dynamics
    s = iso.spin
    P = [];
    norm = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)
    if s.target == "max"
        Pz = -m[4,end]/(s.Nspins*norm)
        P = [0.0, 0.0, 0.0, Pz]
        
    elseif s.target == "min"
        Px = m[2,end]/(s.Nspins*norm)
        Py = m[3,end]/(s.Nspins*norm)
        Pz = m[4,end]/(s.Nspins*norm)
        P = [0.0, Px, Py, Pz]
    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end
    
    return P
end

function grad_saturation_contrast_Mx(iso::Isochromat)
    m = iso.magnetization.dynamics
    s = iso.spin
    P = [];
    norm = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)
    if s.target == "max"
        Px = -m[2,end]/(s.Nspins*norm)
        P = [0.0, Px, 0.0, 0.0]
        
    elseif s.target == "min"
        Px = m[2,end]/(s.Nspins*norm)
        Py = m[3,end]/(s.Nspins*norm)
        Pz = m[4,end]/(s.Nspins*norm)
        P = [0.0, Px, Py, Pz]
    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end
    
    return P
end

function grad_saturation_contrast_square(iso::Isochromat)
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        Pz = -m[4,end]/s.Nspins
        P = [0.0, 0.0, 0.0, Pz]
        
    elseif s.target == "min"
        Px = m[2,end]/s.Nspins
        Py = m[3,end]/s.Nspins
        Pz = m[4,end]/s.Nspins
        P = [0.0, Px, Py, Pz]
    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end
    
    return P
end


function grad_target_steady_state(iso::Isochromat)
    s = iso.spin
    if s.target == "min"
        # Steady State
        ss = steady_state_matrix(s)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

        # Magnetization
        m  = iso.magnetization.dynamics
        Mx = m[2,end]
        My = m[3,end]
        Mz = m[4,end]

        # Adjoint State
        norm = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2)
        Px = (Mx - Mx_ss)/(s.Nspins*norm)
        Py = (My - My_ss)/(s.Nspins*norm)
        Pz = (Mz - Mz_ss)/(s.Nspins*norm)

    elseif s.target == "max"
        # Steady State
        ss = steady_state_matrix(s)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

        # Magnetization
        m  = iso.magnetization.dynamics
        Mx = m[2,end]
        My = m[3,end]
        Mz = m[4,end]

        # Adjoint State
        norm = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)
        Px = (Mx - Mx_ss)/(s.Nspins*norm)
        Py = (My - My_ss)/(s.Nspins*norm)
        Pz = (Mz - Mz_ss)/(s.Nspins*norm)

    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end

    return [0.0, Px, Py, Pz]
end

function grad_targetB0_steady_state(iso::Isochromat)
    s = iso.spin

    # Steady State
    ss = steady_state_matrix(s)
    Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

    # Magnetization
    m  = iso.magnetization.dynamics
    Mx = m[2,end]
    My = m[3,end]
    Mz = m[4,end]

    # Adjoint State
    norm = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)
    Px = (Mx - Mx_ss)/(s.Nspins*norm)
    Py = (My - My_ss)/(s.Nspins*norm)
    Pz = (Mz - Mz_ss)/(s.Nspins*norm)

    return [0.0, Px, Py, Pz]
end



function grad_saturation_contrast_steady_state(iso::Isochromat)
    s = iso.spin

    if s.target == "max"
        # Steady State
        ss = steady_state_matrix(s)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

        # Magnetization
        m  = iso.magnetization.dynamics
        Mx = m[2,end]
        My = m[3,end]
        Mz = m[4,end]

        # Adjoint State
        # norm = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2)
        Px = 0.0 #-(Mx - Mx_ss)/(s.Nspins*norm)
        Py = 0.0 #-(My - My_ss)/(s.Nspins*norm)
        Pz = -(Mz - Mz_ss)/(s.Nspins*sqrt((Mz - Mz_ss)^2 + 1e-15))
    
    elseif s.target == "min"
        # Magnetization
        m  = iso.magnetization.dynamics
        Px = m[2,end]/(s.Nspins*sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15))
        Py = m[3,end]/(s.Nspins*sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15))
        Pz = m[4,end]/(s.Nspins*sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15))

    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end

    return [0.0, Px, Py, Pz]
end