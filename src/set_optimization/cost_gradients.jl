"""
Cost gradients dictonary

"""
#### FUNCTIONS ####                    
function grad_euclidean_norm(iso::Magnetization)
    mag = iso.magnetization
    M = [mag[2,end]; mag[3,end]; mag[4,end]]
    M_norm = norm(M)
    Mx_tf = mag[2,end]/M_norm
    My_tf = mag[3,end]/M_norm
    Mz_tf = mag[4,end]/M_norm

    P = [1.0, Mx_tf, My_tf, Mz_tf]

    return P
end


function grad_target_one_spin(iso::Magnetization; M_tar = [0.5; 0.5; 0.0])    
    # Target Magnetization
    Mx_tar = M_tar[1,1]
    My_tar = M_tar[2,1]
    Mz_tar = M_tar[3,1]

    # Magnetization
    mag = iso.magnetization  
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

function grad_contrast()
end

function grad_euclidean_distance(iso::Magnetization)
    c = 0.0;
        c_max = 0.0;
        c_min = 0.0;
        for i ∈ eachindex(iso.spin)
            mag = iso.magnetization
            spin = iso.spin[i]
            if spin.target == "max"
                c_max = -√(sum(mag[2:end,end].*mag[2:end,end]))
            elseif spin.target == "min"
                c_min = √(sum(mag[2:end,end].*mag[2:end,end]))
            else
                continue
            end
        end
        c = c_max + c_min
    println("Cost: euclidean distance = $(c_max)")
    return c
end

cost_gradients = Dict("Euclidean Norm"  => grad_euclidean_norm,
                      "Target One Spin" => grad_target_one_spin,
                      "Contrast"        => grad_contrast)    

