#### Finite Difference ####
function finite_difference_cost(cost_func::Function, iso::Isochromat, ΔM::Float64)
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
        finite_diffs[i,1] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / ΔM
        
        # Reset the perturbed value for the next iteration
        M_perturbed[i+1,end] = M_perturbed[i+1,end] - ΔM 
    
    end
    return finite_diffs
end


function finite_difference_field(cost_func::Function, cf::ControlField, spin::Spin, Δcf::Float64, field::String)
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
            cf_perturbed = ControlField(perturbation, cf.B1y, cf.B1x_max_amp, 
                            cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)

            # With perturbation
            mag_perturbed = forward_propagation(cf_perturbed, spin)
            dyn_perturbed = Magnetization(mag_perturbed)
            iso_perturbed = Isochromat(dyn_perturbed, spin)

            # Calculate the finite difference for the current variable
            finite_diffs[1, i] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / Δcf

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
            cf_perturbed = ControlField(cf.B1x, perturbation, cf.B1x_max_amp, 
                            cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)

            # With perturbation
            mag_perturbed = forward_propagation(cf_perturbed, spin)
            dyn_perturbed = Magnetization(mag_perturbed)
            iso_perturbed = Isochromat(dyn_perturbed, spin)


            # Calculate the finite difference for the current variable
            finite_diffs[1, i] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / Δcf

            # Reset the perturbed value for the next iteration
            perturbation[1, i] = perturbation[1, i] - Δcf
        end       
    end

    return finite_diffs
end
