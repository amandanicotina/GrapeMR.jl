# Optimizes based on cost function, gradient and propagation
## Optimize function with all inputs for the optimization

function grape_optimize(op::OptimizationParams, cf::ControlFields, spin::Spins; max_iter=2500, ϵ = 1e-4) 
    # Initializing Optimization
    B1x            = copy(cf.B1x);
    B1x_all        = zeros(1, op.N, max_iter+1);
    B1x_all[1,:,1] = B1x;
    B1x_optimize   = op.fields_opt[1]

    B1y            = copy(cf.B1y);
    B1y_all        = zeros(1, op.N, max_iter+1);
    B1y_all[1,:,1] = B1y;
    B1y_optimize   = op.fields_opt[2]

    cost_vals = zeros(1, max_iter+1);
    final_iso = Vector{Magnetization}(undef, 2)

    for k ∈ 1:length(op.fields_opt)
        if B1x_optimize
            println("Entering B1x optimization")
            for i ∈ 1:max_iter
                control_Bx = ControlFields(B1x, B1y, cf.B1x_max_amp,
                                    cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step);
                # Propagation
                mag_Bx = forward_propagation(control_Bx, spin)
                iso_Bx = Magnetization((mag_Bx,), (spin,))

                # Cost function
                cost_func = cost_functions[op.cost_function](iso_Bx)
                cost_vals[1, i] = cost_func
                println("Cost function value = ", cost_func)

                # Control field calculation
                B1x_arr = update_control_field(control_Bx, iso_Bx, Ix, op.cost_function, ϵ)   
                B1x = B1x_arr
                B1x_all[1, :, i+1] = B1x_arr
            end
            B1x_optimize = false

        elseif B1y_optimize
            println("Entering B1y optimization")
            for j ∈ 1:max_iter
                #B1x        = copy(cf.B1x);
                control_By = ControlFields(B1x, B1y, cf.B1x_max_amp,
                                    cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step);
                # Propagation
                mag_By = forward_propagation(control_By, spin)
                iso_By = Magnetization((mag_By,), (spin,))

                # Cost function
                cost_func = cost_functions[op.cost_function](iso_By)
                cost_vals[1, j] = cost_func;
                println("Cost function value = ", cost_func)

                # Control field calculation
                B1y_arr = update_control_field(control_By, iso_By, Iy, op.cost_function, ϵ)   
                B1y = B1y_arr
                B1y_all[1, :, j+1] = B1y_arr
            end
            B1y_optimize = false
        end
    end

    # Calculating with optimized fields
    all_controls  = ControlFields(B1x_all, B1y_all, cf.B1x_max_amp, cf.B1y_max_amp, 
                                    cf.t_control, cf.band_width, cf.band_width_step)
    final_control = ControlFields(B1x_all[:,:,end], B1y_all[:,:,end], cf.B1x_max_amp, cf.B1y_max_amp, 
                                    cf.t_control, cf.band_width, cf.band_width_step)

    mag_opt       = forward_propagation(final_control, spin)
    final_iso     = Magnetization((mag_opt,), (spin,))

    return final_iso, all_controls
end 




#### Finite Difference ####
function finite_difference_cost(cost_func::Function, iso::Magnetization, ΔM::Float64)
    # Make a copy of the variable values to avoid modifying the original
    spin        = iso.spin
    iso_vals    = iso
    M_perturbed = copy(iso.magnetization[1][:,end])

    # Initialize an array to store finite differences for each variable
    finite_diffs = zeros(Float64, 3, 1)
    for i ∈ 1:3
        M_perturbed[i+1,end] = M_perturbed[i+1,end] + ΔM
        iso_perturbed = Magnetization((M_perturbed,), spin,)

        # Calculate the finite difference for the current variable
        finite_diffs[i,1] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / ΔM
        
        # Reset the perturbed value for the next iteration
        M_perturbed[i+1,end] = M_perturbed[i+1,end] - ΔM 
    
    end
    return finite_diffs
end


function finite_difference_field(cost_func::Function, cf::ControlFields, spin::Spins, Δcf::Float64)
    finite_diffs = zeros(Float64, 1, length(cf.B1x))

    # Copy of the variable values to avoid modifying the original
    cf_vals      = cf.B1x
    perturbation = copy(cf_vals) 

    # No perturbation
    mag_vals = forward_propagation(cf, spin)
    iso_vals = Magnetization((mag_vals,), (spin,))

    for i ∈ 1:length(cf.B1x)
        perturbation[1, i] = perturbation[1, i] + Δcf
        cf_perturbed = ControlFields(cf.B1x, perturbation, cf.B1x_max_amp, 
                        cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)

        # With perturbation
        mag_perturbed = forward_propagation(cf_perturbed, spin)
        iso_perturbed = Magnetization((mag_perturbed,), (spin,))

        # Calculate the finite difference for the current variable
        finite_diffs[1, i] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / Δcf

        # Reset the perturbed value for the next iteration
        perturbation[1, i] = perturbation[1, i] - Δcf
    end

    return finite_diffs
end
