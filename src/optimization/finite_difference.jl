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
        (cost_pert, grad) = GrapeMR.cost_function(iso_perturbed, cost)
        (cost_vals, grad) = GrapeMR.cost_function(iso_vals, cost)
        finite_diffs[i,1] = (cost_pert - cost_vals) / ΔM
        
        # Reset the perturbed value for the next iteration
        M_perturbed[i+1,end] = M_perturbed[i+1,end] - ΔM 
    
    end
    return finite_diffs
end

function finite_difference_field(cost::Symbol, cf::ControlField, spin::Spins, Δcf::Float64, field::String)
    finite_diffs = zeros(Float64, 1, length(cf.B1x))

    # Get the original isochromat values (no perturbation)
    iso_vals = dynamics(cf, spin)

    # Helper function for applying perturbations and calculating finite differences
    function calculate_finite_difference(cf::ControlField, spin::Spins, perturbation, index::Int, Δcf::Float64, field::String)
        if field == "B1x"
            # Perturb B1x field
            perturbed_cf_pos = ControlField(copy(perturbation), cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)
            perturbation[1, index] -= 2 * Δcf  
            perturbed_cf_neg = ControlField(copy(perturbation), cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)
        elseif field == "B1y"
            # Perturb B1y field
            perturbed_cf_pos = ControlField(cf.B1x, copy(perturbation), cf.B1_ref, cf.Bz, cf.t_control)
            perturbation[1, index] -= 2 * Δcf  # Perturb in the negative direction
            perturbed_cf_neg = ControlField(cf.B1x, copy(perturbation), cf.B1_ref, cf.Bz, cf.t_control)
        else
            error("Parameter not defined. Acceptable inputs are \"B1x\" or \"B1y\"")
        end
        
        # Dynamics with positive perturbation
        iso_pos = dynamics(perturbed_cf_pos, spin)
        
        # Dynamics with negative perturbation
        iso_neg = dynamics(perturbed_cf_neg, spin)
        
        # Central difference formula
        return (cost_function(iso_pos, cost)[1] - cost_function(iso_neg, cost)[1]) / (2 * Δcf)
    end

    # Select the appropriate field to perturb
    if field == "B1x"
        perturbation = copy(cf.B1x)
        for i ∈ 1:length(cf.B1x)
            perturbation[1, i] += Δcf  # Perturb in the positive direction
            finite_diffs[1, i] = calculate_finite_difference(cf, spin, perturbation, i, Δcf, "B1x")
            perturbation[1, i] -= Δcf  # Reset perturbation
        end
    elseif field == "B1y"
        perturbation = copy(cf.B1y)
        for i ∈ 1:length(cf.B1y)
            perturbation[1, i] += Δcf  # Perturb in the positive direction
            finite_diffs[1, i] = calculate_finite_difference(cf, spin, perturbation, i, Δcf, "B1y")
            perturbation[1, i] -= Δcf  # Reset perturbation
        end
    else
        error("Parameter not defined. Acceptable inputs are \"B1x\" or \"B1y\"")
    end

    return finite_diffs
end
