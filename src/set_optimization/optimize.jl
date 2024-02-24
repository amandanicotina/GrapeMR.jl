# Optimizes based on cost function, gradient and propagation
## Optimize function with all inputs for the optimization

struct GrapeOutput
    magnetization::Vector{Magnetization}
    control_field::ControlFields
end


function grape_optimize(op::OptimizationParams, cf::ControlFields, spins::Array{Spins}; max_iter=2500, ϵ = 1e-4) 
    # Initializing Optimization
    u1x            = copy(cf.B1x);
    u1x_all        = zeros(1, op.N, max_iter+1);
    u1x_all[1,:,1] = u1x;

    u1y            = copy(cf.B1y);
    u1y_all        = zeros(1, op.N, max_iter+1);
    u1y_all[1,:,1] = u1y;

    # TODO: map instead
    for spin ∈ spins
        update_control_fields(spin, op, cf, max_iter, u1x, u1x_all, u1y, u1y_all, ϵ)
    end

    # Calculating with optimized fields
    all_controls  = ControlFields(u1x_all, u1y_all, cf.B1x_max_amp, cf.B1y_max_amp, 
                                    cf.t_control, cf.band_width, cf.band_width_step)

    final_control = ControlFields(u1x_all[:,:,end], u1y_all[:,:,end], cf.B1x_max_amp, cf.B1y_max_amp, 
                                    cf.t_control, cf.band_width, cf.band_width_step)

    function final_iso(spin::Spins)
        mag_opt       = forward_propagation(final_control, spin)
        return Magnetization(mag_opt, spin)
    end

    return GrapeOutput(map(final_iso, spins), all_controls)
end 

function update_control_fields(spin::Spins, params::OptimizationParams, cf::ControlFields, max_iter::Int64, u1x::AbstractArray, u1x_all::Array{Float64}, u1y::AbstractArray, u1y_all::Array{Float64}, ϵ::Float64)
    cost_func = cost_functions[params.cost_function]
    cost_vals = zeros(1, max_iter+1);
    if params.fields_opt[1]
        println("Entering B1x optimization")
        for i ∈ 1:max_iter
            u1y        = copy(cf.B1y);
            control_ux = ControlFields(u1x, u1y, cf.B1x_max_amp,
                                cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step);
            # Propagation
            mag_ux = forward_propagation(control_ux, spin)
            iso_ux = Magnetization(mag_ux, spin)

            # Cost function
            cost_vals[1, i] = cost_func(iso_ux)
            println("Cost function value = ", cost_vals[1, i])

            # Control field calculation
            u1x_arr = update_control_field(control_ux, iso_ux, Ix, params.cost_function, ϵ)   
            u1x = u1x_arr
            u1x_all[1, :, i+1] = u1x_arr
        end
    end

    if params.fields_opt[2]
        println("Entering B1y optimization")
        for j ∈ 1:max_iter
            u1x        = copy(cf.B1x);
            control_uy = ControlFields(u1x, u1y, cf.B1x_max_amp,
                                cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step);
            # Propagation
            mag_uy = forward_propagation(control_uy, spin)
            iso_uy = Magnetization(mag_uy, spin)

            # Cost function
            cost_vals[1, j] = cost_func(iso_uy);
            println("Cost function values = ", cost_vals[1, j])

            # Control field calculation
            u1y_arr = update_control_field(control_uy, iso_uy, Iy, params.cost_function, ϵ)   
            u1y = u1y_arr
            u1y_all[1, :, j+1] = u1y_arr
        end
    end
end

#### Finite Difference ####
function finite_difference_cost(cost_func::Function, iso::Magnetization, ΔM::Float64)
    # Make a copy of the variable values to avoid modifying the original
    spin        = iso.spin
    iso_vals    = iso
    M_perturbed = copy(iso.magnetization[:,end])

    # Initialize an array to store finite differences for each variable
    finite_diffs = zeros(Float64, 3, 1)
    for i ∈ 1:3
        M_perturbed[i+1,end] = M_perturbed[i+1,end] + ΔM
        iso_perturbed = Magnetization(M_perturbed, spin)

        # Calculate the finite difference for the current variable
        finite_diffs[i,1] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / ΔM
        
        # Reset the perturbed value for the next iteration
        M_perturbed[i+1,end] = M_perturbed[i+1,end] - ΔM 
    
    end
    return finite_diffs
end


function finite_difference_field(cost_func::Function, cf::ControlFields, spin::Spins, Δcf::Float64, field::String)
    finite_diffs = zeros(Float64, 1, length(cf.B1x))

    if field == "B1x"
        # Copy of the variable values to avoid modifying the original
        cf_vals      = cf.B1x
        perturbation = copy(cf_vals) 

        # No perturbation
        mag_vals = forward_propagation(cf, spin)
        iso_vals = Magnetization(mag_vals, spin)

        for i ∈ 1:length(cf.B1x)
            perturbation[1, i] = perturbation[1, i] + Δcf
            cf_perturbed = ControlFields(perturbation, cf.B1y, cf.B1x_max_amp, 
                            cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)

            # With perturbation
            mag_perturbed = forward_propagation(cf_perturbed, spin)
            iso_perturbed = Magnetization(mag_perturbed, spin)

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
        iso_vals = Magnetization(mag_vals, spin)

        for i ∈ 1:length(cf.B1y)
            perturbation[1, i] = perturbation[1, i] + Δcf
            cf_perturbed = ControlFields(cf.B1x, perturbation, cf.B1x_max_amp, 
                            cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)

            # With perturbation
            mag_perturbed = forward_propagation(cf_perturbed, spin)
            iso_perturbed = Magnetization(mag_perturbed, spin)

            # Calculate the finite difference for the current variable
            finite_diffs[1, i] = (cost_func(iso_perturbed) - cost_func(iso_vals)) / Δcf

            # Reset the perturbed value for the next iteration
            perturbation[1, i] = perturbation[1, i] - Δcf
        end       
    end

    return finite_diffs
end
