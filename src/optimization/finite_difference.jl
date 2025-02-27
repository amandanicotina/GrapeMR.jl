#### Finite Difference ####
function finite_difference_cost(cost::Function, iso::Isochromat, ΔM::Float64)
    # Make a copy of the variable values to avoid modifying the original
    spin = iso.spin
    iso_vals = iso
    M_perturbed = copy(iso.magnetization.dynamics[:, end])

    # Initialize an array to store finite differences for each variable
    finite_diffs = zeros(Float64, 3, 1)
    for i ∈ 1:3
        M_perturbed[i+1, end] = M_perturbed[i+1, end] + ΔM
        dyn_perturbed = Magnetization(M_perturbed)
        iso_perturbed = Isochromat(dyn_perturbed, spin)

        # Calculate the finite difference for the current variable
        (cost_pert, grad) = cost(iso_perturbed)
        (cost_vals, grad) = cost(iso_vals)
        finite_diffs[i, 1] = (cost_pert - cost_vals) / ΔM

        # Reset the perturbed value for the next iteration
        M_perturbed[i+1, end] = M_perturbed[i+1, end] - ΔM

    end
    return finite_diffs
end


function finite_difference_field(s::Spin, cf::ControlField, gp::GrapeParams, field::String, Δcf::Float64)
    finite_difference = zeros(eltype(cf.B1x), 1, length(cf.B1x))
    for i ∈ 1:length(cf.B1x)
        perturbation_pos = copy(cf.B1x)
        perturbation_neg = copy(cf.B1x)
        perturbation_pos[1, i] += Δcf
        perturbation_neg[1, i] -= Δcf
        if field == "B1x"
            perturbed_cf_pos = ControlField(copy(perturbation_pos), cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)
            perturbed_cf_neg = ControlField(copy(perturbation_neg), cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)
        elseif field == "B1y"   
            perturbed_cf_pos = ControlField(cf.B1x, copy(perturbation_pos), cf.B1_ref, cf.Bz, cf.t_control)
            perturbed_cf_neg = ControlField(cf.B1x, copy(perturbation_neg), cf.B1_ref, cf.Bz, cf.t_control)
        else
            error("Parameter not defined. Acceptable inputs are \"B1x\" or \"B1y\"")
        end
        iso_pos = dynamics(perturbed_cf_pos, s) 
        iso_neg = dynamics(perturbed_cf_neg, s)
        cost_vars_pos = gp.cost_function(iso_pos)
        cost_pos, _ = cost_vars_pos
        cost_vars_neg = gp.cost_function(iso_neg)
        cost_neg, _ = cost_vars_neg   
        finite_difference[1, i] = (cost_pos - cost_neg) / (2 * Δcf)

        # Reset perturbation
        perturbation_pos[1, i] -= Δcf
        perturbation_neg[1, i] += Δcf
    end
    return finite_difference
end
    