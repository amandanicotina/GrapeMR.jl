# Optimizes based on cost function, gradient and propagation
## Optimize function with all inputs for the optimization

struct GrapeOutput
    isochromats::Vector{Isochromat}
    control_field::ControlField
    cost_values::Array{Float64}
end


function grape_optimize(op::OptimizationParams, cf::ControlField, spins::Vector{Spin}; max_iter=2500, ϵ = 1e-4)
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
    all_controls  = ControlField(u1x_all, u1y_all, cf.B1x_max_amp, cf.B1y_max_amp, 
                                    cf.t_control, cf.band_width, cf.band_width_step)

    final_control = ControlField(u1x_all[:,:,end], u1y_all[:,:,end], cf.B1x_max_amp, cf.B1y_max_amp, 
                                    cf.t_control, cf.band_width, cf.band_width_step)

    function final_iso(spin::Spin)
        mag_opt       = forward_propagation(final_control, spin)
        return Magnetization(mag_opt)
    end

    return GrapeOutput(map(final_iso, spins), all_controls)
end 

function update_control_fields(spin::Spin, params::OptimizationParams, cf::ControlField, max_iter::Int64, u1x::AbstractArray, u1x_all::Vector{Float64}, u1y::AbstractArray, u1y_all::Vector{Float64}, ϵ::Float64)
    cost_func = cost_functions[params.cost_function]
    cost_vals = zeros(1, max_iter+1);
    if params.fields_opt[1]
        println("Entering B1x optimization")
        for i ∈ 1:max_iter
            u1y        = copy(cf.B1y);
            control_ux = ControlField(u1x, u1y, cf.B1x_max_amp,
                                cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step);
            # Propagation
            mag_ux = forward_propagation(control_ux, spin)
            iso_ux = Magnetization(mag_ux)

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
            control_uy = ControlField(u1x, u1y, cf.B1x_max_amp,
                                cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step);
            # Propagation
            mag_uy = forward_propagation(control_uy, spin)
            iso_uy = Magnetization(mag_uy)

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



function grape(op::OptimizationParams, cf::ControlField, spins::Vector{Spin}; max_iter=2500, ϵ = 1e-4)
    b1x_old = copy(cf.B1x)
    b1y_old = copy(cf.B1y)
    ∇x = zeros(Float64, 1, op.N)
    ∇y = zeros(Float64, 1, op.N)
    #cost_vals = zeros(Float64, length(spins), max_iter)
    cost_vals = zeros(Float64, 2, length(spins), max_iter)
    grape_output = GrapeOutput([], cf, cost_vals)

    for i ∈ 1:max_iter

        for (j, spin) ∈ enumerate(spins)
            # Propagation
            mag = forward_propagation(cf, spin)
            dyn = Magnetization(mag)
            iso = Isochromat(dyn, spin)
            adj = backward_propagation(cf, iso, grad_euclidean_norm)
            push!(grape_output.isochromats, iso)
            # Gradient
            if op.fields_opt[1]
                ∇x += gradient(adj, mag, Ix)
                (xu1x, _) = update(cf, (∇x, zeros(Float64, 1, op.N)), ϵ)
                cfx = ControlField(xu1x, b1y_old, cf.B1x_max_amp, cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)
                magx = forward_propagation(cfx, spin)
                dynx = Magnetization(magx)
                isox = Isochromat(dynx, spin)
                grape_output.cost_values[1,j,i] = op.cost_function(isox)
            end
            if op.fields_opt[2]
                ∇y += gradient(adj, mag, Iy)
                (_, yu1y) = update(cf, (zeros(Float64, 1, op.N), ∇y), ϵ)
                cfy = ControlField(b1x_old, yu1y, cf.B1x_max_amp, cf.B1y_max_amp, cf.t_control, cf.band_width, cf.band_width_step)
                magy = forward_propagation(cfy, spin)
                dyny = Magnetization(magy)
                isoy = Isochromat(dyny, spin)
                grape_output.cost_values[2,j,i] = op.cost_function(isoy)
            end 
        end

        # Control Field
        (u1x, u1y) = update(cf, (∇x, ∇y), ϵ)
        cf.B1x = u1x
        cf.B1y = u1y
    end

    return grape_output
end


"""
    gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)

gradient
    # Input  
    - χ = (::Matrix{Float64}) - Adjoint State
    - M = (::Matrix{Float64}) - Forward Propagation
    - H = (::Matrix) - Hamiltonian

    # Output
    - ΔJ - 1xN matrix
"""
function gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)
    grad = zeros(Float64, 1, length(M[1,:])-1)
    for i ∈ 1:(length(M[1,:])-1)
        grad[1,i] = transpose(χ[:,i+1])*H*M[:,i]./2π
    end
    return grad
end



"""
    update(cf::InitialControlFields, ∇xy::gradient, ϵ::Float64)

update
    # Input  
    - cf:  (::InitialControlFields) - Control fields struct
    - iso: (::Magnetization) -
    - cost_function: (::String) - Key to the cost function dictionary
    - ϵ:   (::Float64) - 

    # Output
    - Control Field - 1xN matrix
"""
function update(cf::ControlField, ∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ*∇xy[1]
    u1y = cf.B1y .- ϵ*∇xy[2]
    return u1x, u1y
end



