#### Finite Difference ####
function finite_difference_cost(cost::Symbol, iso::Isochromat, ΔM::Float64)
    # Make a copy of the variable values to avoid modifying the original
    spin        = iso.spin
    iso_vals    = iso
    M_perturbed = copy(iso.magnetization.dynamics[:,end])

    # Initialize an array to store finite differences for each variable
    finite_diffs = zeros(Float64, 3, 1)
    for i ∈ 1:3
        M_perturbed[i+1,end] = M_perturbed[i+1,end] + ΔM
        dyn_perturbed = Magnetization(M_perturbed)
        iso_perturbed = Isochromat(dyn_perturbed, spin)

        # Calculate the finite difference for the current variable
        finite_diffs[i,1] = (GrapeMR.cost_function(iso_perturbed, cost) - GrapeMR.cost_function(iso_vals, cost)) / ΔM
        
        # Reset the perturbed value for the next iteration
        M_perturbed[i+1,end] = M_perturbed[i+1,end] - ΔM 
    
    end
    return finite_diffs
end


function finite_difference_field(cost::Symbol, cf::ControlField, spin::Spins, Δcf::Float64, field::String)
    finite_diffs = zeros(Float64, 1, length(cf.B1x))

    if field == "B1x"
        # Copy of the variable values to avoid modifying the original
        cf_vals      = cf.B1x
        perturbation = copy(cf_vals) 

        # No perturbation
        mag_vals = forward_propagation(cf, spin)
        dyn_vals = Magnetization(mag_vals)
        iso_vals = Isochromat(dyn_vals, spin)

        for i ∈ 1:length(cf.B1x)
            perturbation[1, i] = perturbation[1, i] + Δcf
            cf_perturbed = ControlField(perturbation, cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)

            # With perturbation
            mag_perturbed = forward_propagation(cf_perturbed, spin)
            dyn_perturbed = Magnetization(mag_perturbed)
            iso_perturbed = Isochromat(dyn_perturbed, spin)

            # Calculate the finite difference for the current variable
            finite_diffs[1, i] = (cost_function(iso_perturbed, cost) - cost_function(iso_vals, cost)) / Δcf

            # Reset the perturbed value for the next iteration
            perturbation[1, i] = perturbation[1, i] - Δcf
        end

    else
        # Copy of the variable values to avoid modifying the original
        cf_vals      = cf.B1y
        perturbation = copy(cf_vals) 

        # No perturbation
        mag_vals = forward_propagation(cf, spin)
        dyn_vals = Magnetization(mag_vals)
        iso_vals = Isochromat(dyn_vals, spin)

        for i ∈ 1:length(cf.B1y)
            perturbation[1, i] = perturbation[1, i] + Δcf
            cf_perturbed = ControlField(cf.B1x, perturbation, cf.B1_ref, cf.Bz, cf.t_control)

            # With perturbation
            mag_perturbed = forward_propagation(cf_perturbed, spin)
            dyn_perturbed = Magnetization(mag_perturbed)
            iso_perturbed = Isochromat(dyn_perturbed, spin)


            # Calculate the finite difference for the current variable
            finite_diffs[1, i] = (cost_function(iso_perturbed, cost) - cost_function(iso_vals, cost)) / Δcf

            # Reset the perturbed value for the next iteration
            perturbation[1, i] = perturbation[1, i] - Δcf
        end       
    end

    return finite_diffs
end
