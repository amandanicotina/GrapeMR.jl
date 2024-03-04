"""
Cost Functions dictonary

long  = range(-1, 1, length=N);
trans = range(-1, 1, length=N);

cost_func(Mx, My, Mz) = sqrt.(Mz.^2 .+ Mt.^2);
# Plot cost function values
contourf(trans, long, cost_func, color = :jet)

"""


#### FUNCTIONS ####                    
function cost_euclidean_norm(iso::Isochromat)
    mag = iso.magnetization.dynamics
    Mx = mag[2,end]
    My = mag[3,end]
    Mz = mag[4,end]

    J = (Mx^2 + My^2 + Mz^2)/2

    return J
end

function cost_target_one_spin(iso::Isochromat; M_tar = [0.5; 0.5; 0.0])
    J_tar = 0
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


function cost_saturation_contrast(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin
    if s.target == "max"
        c = 1/2 - sqrt(m[4,end]*m[4,end] + 1e-15)/2
    elseif s.target == "min"
        c = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)/2
    end
    return c
end

function cost_saturation_contrast_square(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin
    if s.target == "max"
        c = 1/2 - sum(m[4,end]*m[4,end])/2
        # c = - sum(m[4,end]*m[4,end])/2
    elseif s.target == "min"
        c = sum(m[2:end,end].*m[2:end,end])/2
    end
    return c
end

function cost_steady_state(iso::Isochromat, ss::Tuple)
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


#cost_functions = Dict("Euclidean Norm"      => euclidean_norm,
#                      "Target One Spin"     => target_one_spin,
#                      "Saturation Contrast" => saturation_contrast)    

