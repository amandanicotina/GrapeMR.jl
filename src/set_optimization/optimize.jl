# Optimizes based on cost function, gradient and propagation
## Optimize function with all inputs for the optimization

function grape_optimize(cf::InitialControlFields, spin::Spins; max_iter=1500, ϵ = 1e-7) 
    B1 = cf.B1x_init_control;
    B1_all = zeros(1, cf.N, max_iter+1);
    B1_all[1,:,1] = B1;

    cost_vals = zeros(1, max_iter+1);
    isoList = [];
    for i ∈ 1:max_iter
        init_control_field = InitialControlFields(cf.N, B1, cf.B1x_max_amplitude, cf.B1y_init_control,
                            cf.B1y_max_amplitude, cf.t_control, 0.0, 0.0);
        # Propagation
        mag = magnetization_ODE(init_control_field, spin)
        isoList = Magnetization((mag,), (spin,))

        # Cost function
        cost_func = cost_functions["Euclidean Norm"](isoList)
        cost_vals[1, i] = cost_func
        println("Cost function value = ", cost_func)

        # Control field calculation
        B1_arr_op = update_control_field(init_control_field, spin, isoList, ϵ)  
        B1 = B1_arr_op
        B1_all[1, :, i+1] = B1_arr_op

    end

   # return in a struct -> controls = ControlFields(B1_all, ...)
    return isoList, B1_all, cost_vals
end




#### Finite Difference ####
function finite_difference_cost(cost_func::Function, iso::Magnetization, ΔM::Float64)
    # Make a copy of the variable values to avoid modifying the original
    spin = iso.spin
    iso_vals = iso
    M_perturbed = iso.magnetization[1][:,end]

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


function finite_difference_field(cost_func::Function, cf::InitialControlFields, spin::Spins, Δcf::Float64)
    finite_diffs = zeros(Float64, 1, cf.N)

    # Copy of the variable values to avoid modifying the original
    cf_vals      = cf.B1x_init_control
    perturbation = copy(cf_vals) 

    # No perturbation
    mag_vals = magnetization_ODE(cf, spin)
    iso_vals = Magnetization((mag_vals,), (spin,))

    for i ∈ 1:cf.N
        perturbation[1, i] = perturbation[1, i] + Δcf
        cf_perturbed = InitialControlFields(cf.N, perturbation, cf.B1x_max_amplitude, cf.B1y_init_control,
                        cf.B1y_max_amplitude, cf.t_control, cf.band_width, cf.band_width_step)

        # With perturbation
        mag_perturbed = magnetization_ODE(cf_perturbed, spin)
        iso_perturbed = Magnetization((mag_perturbed,), (spin,))

        # Calculate the finite difference for the current variable
        finite_diffs[1, i] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / Δcf

       # println("Finite Diff val = ", finite_diffs)
        #println("Cost vals = ", cost_func(iso_vals))
        #println("Cost perturbed = ", cost_func(iso_perturbed))

        # Reset the perturbed value for the next iteration
        perturbation[1, i] = perturbation[1, i] - Δcf
    end

    return finite_diffs
end
