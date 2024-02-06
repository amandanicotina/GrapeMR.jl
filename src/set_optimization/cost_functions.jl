"""
Cost Functions dictonary

long  = range(-1, 1, length=N);
trans = range(-1, 1, length=N);

cost_func(Mx, My, Mz) = sqrt.(Mz.^2 .+ Mt.^2);
# Plot cost function values
contourf(trans, long, cost_func, color = :jet)

"""
#### FUNCTIONS ####                    
function euclidean_norm(iso::Magnetization)
    mag = iso.magnetization[1]
    Mx = mag[2,end]
    My = mag[3,end]
    Mz = mag[4,end]

    J = sqrt(Mx^2 + My^2 + Mz^2)

    return J
end

function target_one_spin()
    
end

function contrast()
    
end

function euclidean_distance(iso::Magnetization)
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

cost_functions = Dict("Euclidean Norm"  => euclidean_norm,
                      "Target One Spin" => target_one_spin,
                      "Contrast"        => contrast)    

