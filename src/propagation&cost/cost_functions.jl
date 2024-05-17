"""
Cost Functions 

"""
function cost_function(iso::Isochromat, cf::Symbol)
    @match cf begin
        :euclidean_norm      => euclidean_norm(iso)
        :target_one_spin     => target_one_spin(iso)
        :target_steady_state => target_steady_state(iso, ())
        #:target_steady_state_opt    => target_steady_statee_opt()
        :saturation_contrast => saturation_contrast(iso)
        :saturation_contrast_square => saturation_contrast_square(iso)
        _                    => error("Cost function not defined")
    end
end

              
function euclidean_norm(iso::Isochromat)
    mag = iso.magnetization.dynamics
    Mx = mag[2,end]
    My = mag[3,end]
    Mz = mag[4,end]

    J = (Mx^2 + My^2 + Mz^2)/2

    return J
end

function target_one_spin(iso::Isochromat; M_tar = [0.0; 1.0; 0.0])

    # Target Magnetization
    Mx_tar = M_tar[1,1]
    My_tar = M_tar[2,1]
    Mz_tar = M_tar[3,1]

    # Magnetization
    mag = iso.magnetization.dynamics
    Mx = mag[2,end]
    My = mag[3,end]
    Mz = mag[4,end]

    J_tar = (Mx - Mx_tar)^2 + (My - My_tar)^2 + (Mz - Mz_tar)^2

    return J_tar
end

function saturation_contrast(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        c = (1 - sqrt(m[4,end]*m[4,end] + 1e-15))/s.Nspins
    elseif s.target == "min"
        c = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)/s.Nspins
    end
    return c
end

function saturation_contrast_square(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        c = (1 - sum(m[4,end]*m[4,end]))/s.Nspins
    elseif s.target == "min"
        c = sum(m[2:end,end].*m[2:end,end])/s.Nspins
    end
    return c
end

function target_steady_state(iso::Isochromat, ss::Tuple)
    c_ss = 0
    # Steady State
    Mxy_ss = ss[1];
    Mz_ss  = ss[2];

    # Magnetization
    mag = iso.magnetization.dynamics
    Mxy = (mag[2,end]).^2 .+ (mag[3,end]).^2
    Mz  = mag[4,end]

    c_ss = sqrt((Mxy - Mxy_ss)^2 + (Mz - Mz_ss)^2)
    return c_ss
end

  
function target_steady_state_opt()
    
end
