"""
Cost gradients dictonary

"""
#### FUNCTIONS ####                    
function grad_euclidean_norm(iso::Magnetization)
    mag = iso.magnetization[1]
    M = [mag[2,end]; mag[3,end]; mag[4,end]]
    M_norm = norm(M)
    Mx_tf = mag[2,end]/M_norm
    My_tf = mag[3,end]/M_norm
    Mz_tf = mag[4,end]/M_norm

    P = [0.0, Mx_tf, My_tf, Mz_tf]

    return P
end


function grad_target_one_spin()    
end

function grad_contrast()
end

function grad_euclidean_distance(iso::Magnetization)
    c = 0.0;
        c_max = 0.0;
        c_min = 0.0;
        for i ∈ eachindex(iso.spin)
            mag = iso.magnetization[i]
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

cost_gradients = Dict("Grad Euclidean Norm"  => grad_euclidean_norm,
                      "Grad Target One Spin" => grad_target_one_spin,
                      "Grad Contrast"        => grad_contrast)    

