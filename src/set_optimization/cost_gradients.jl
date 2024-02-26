"""
Cost gradients dictonary

"""
#### FUNCTIONS ####                    
function grad_euclidean_norm(iso::Isochromat)
    mag = iso.magnetization.dynamics
    M = [mag[2,end]; mag[3,end]; mag[4,end]]
    M_norm = norm(M)
    Mx_tf = mag[2,end]
    My_tf = mag[3,end]
    Mz_tf = mag[4,end]

    P = [0.0, Mx_tf, My_tf, Mz_tf]

    return P
end


function grad_target_one_spin(iso::Isochromat; M_tar = [0.5; 0.5; 0.0])    
    # Target Magnetization
    Mx_tar = M_tar[1,1]
    My_tar = M_tar[2,1]
    Mz_tar = M_tar[3,1]

    # Magnetization
    mag = iso.magnetization.dynamics
    Mx  = mag[2,end]
    My  = mag[3,end]
    Mz  = mag[4,end]
    M_norm = sqrt((Mx - Mx_tar)^2 + (My - My_tar)^2 + (Mz - Mz_tar)^2) 

    Mx_tf = (Mx - Mx_tar)/M_norm
    My_tf = (My - My_tar)/M_norm
    Mz_tf = (Mz - Mz_tar)/M_norm

    P_tar = [1.0, Mx_tf, My_tf, Mz_tf]
    return P_tar
end


function grad_saturation_contrast(iso::Isochromat)
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        Pz = -m[4,end]
        P = [0.0, 0.0, 0.0, Pz]
        
    elseif s.target == "min"
        Px = m[2,end]
        Py = m[3,end]
        Pz = m[4,end]
        P = [0.0, Px, Py, Pz]
    end
    
    return P
end

# cost_gradients = Dict("Euclidean Norm"      => grad_euclidean_norm,
#                       "Target One Spin"     => grad_target_one_spin,
#                       "Saturation Contrast" => grad_saturation_contrast)    

