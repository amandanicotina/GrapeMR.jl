const wandb_project::String = "GrapeMR"
Base.broadcastable(cf::ControlField) = Ref(cf)


struct GrapeOutput
    isochromats::Vector{Isochromat}
    control_field::ControlField
    cost_values::Union{Array{Float64}, Metal.MtlVector{Float32, Metal.SharedStorage}}
    params::Parameters
end


"""
    grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})

Implements the grape algorithm [khaneja2005optimal](@cite) 

# Arguments
- `op::OptimizationParams`: Parameters for the optimization itself: max iterations, 
- `gp::GrapeParams`: Parameters related to Grape itself: time points, cost function, mask for which fields are being optimized.
- `cf::ControlField`: Initial control field - spline function -
- `spins::Vector{<:Spins}`: Vector with all spins included in the optimization

# Outputs
A scruct cointaing all optimization results:
- `grape_output::GrapeOutput': Data type with the optimized control fields, spin information and spin dynamics.
"""
function grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter+1) 
    cost_vals    = zeros(Float64, op.max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
    
    for (ϵ, i) ∈ zip(lr_scheduler, 1:op.max_iter)
        # ϵ   = max(ϵ, eps)
        ∇x  = zeros(Float64, 1, gp.N)
        ∇y  = zeros(Float64, 1, gp.N)
        
        for spin ∈ spins
            # Forward Propagation 
            iso = dynamics(cf, spin)
            mag = iso.magnetization.dynamics
            cost_vars = GrapeMR.cost_function(iso, gp.cost_function)
            # Cost Variables
            cost = getindex(cost_vars, 1)
            adj_ini = getindex(cost_vars, 2)
            grape_output.cost_values[i,1] += cost
            # Adjoint Propagation
            adj = backward_propagation(cf, iso, adj_ini)
            
            # Save Isochromats from the last iterations
            if i == op.max_iter
                push!(grape_output.isochromats, iso)
            end
            
            # Gradient
            if gp.fields_opt[1]
                ∇x += gradient(adj, mag, Ix)
            end 
            if gp.fields_opt[2]
                ∇y += gradient(adj, mag, Iy)
            end 
        end

        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x = u1x
        cf.B1y = u1y
    end

    grape_output.control_field.B1x = u1x
    grape_output.control_field.B1y = u1y

    # Print Infos
    final_cost = round(grape_output.cost_values[end], digits = 3)
    println("\n Final Cost Function Value = $final_cost \n")
    RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


function metal_grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter+1) 
    cost_vals    = Metal.zeros(Float32, op.max_iter, 1; storage=Metal.SharedStorage)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
    ∇x  = Metal.zeros(Float32, 1, gp.N; storage=Metal.SharedStorage)
    ∇y  = Metal.zeros(Float32, 1, gp.N; storage=Metal.SharedStorage)
    tmp_costs = Metal.zeros(Float32, length(spins))

    isos = Vector{Isochromat}(undef, length(spins))
    adjs = Vector{Metal.Matrix{Float32}}(undef, length(spins))

    function process_spin(spin_idx::Int64, spin::Spin)
        isos[spin_idx]      = dynamics(cf, spin)
        cost_vars    = GrapeMR.cost_function(isos[spin_idx], gp.cost_function)
        adj_ini      = last(cost_vars)    
        adjs[spin_idx]      = backward_propagation(cf, isos[spin_idx], adj_ini)          
        tmp_costs[spin_idx] = first(cost_vars)
    end
    
    for (ϵ, i) ∈ zip(lr_scheduler, 1:op.max_iter)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)
        fill!(tmp_costs, 0.0)
        # Propagation & cost
        
        @metal threads=length(spins) foreach(process_spin, 1:length(spins), spins)  

        grape_output.cost_values[i, 1] = sum(tmp_costs)
        # Save final magnetization trajectory
        if i == max_iter
            append!(grape_output.isochromats, isos)
        end
        
        # Gradients
        mags = getfield.(isos, :magnetization)
        dyns = getfield.(mags, :dynamics)
        if gp.fields_opt[1]
            ∇x = sum(gradient.(adjs, dyns, Ref(Ix)))
        end 
        if gp.fields_opt[2]
            ∇y = sum(metal_gradient.(adjs, dyns, Ref(Iy)))
        end 
        
        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x = u1x
        cf.B1y = u1y
    end

    grape_output.control_field.B1x = u1x
    grape_output.control_field.B1y = u1y

    # Print Infos
    final_cost = round(grape_output.cost_values[end], digits = 3)
    println("\n Final Cost Function Value = $final_cost \n")
    RF_pulse_analysis(grape_output.control_field)

    return grape_output
end

function threads_grape_metal(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    max_iter     = op.max_iter
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=max_iter + 1) 
    cost_vals    = zeros(Float64, max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)

    # Prepare Metal buffers
    ∇x_gpu = Metal.Buffer(zeros(Float32, 1, gp.N))
    ∇y_gpu = Metal.Buffer(zeros(Float32, 1, gp.N))
    cost_gpu = Metal.Buffer(zeros(Float32, length(spins)))
    isos_gpu = Metal.Buffer(Vector{Isochromat}(undef, length(spins)))
    adjs_gpu = Metal.Buffer(Vector{Matrix{Float32}}(undef, length(spins)))

    for (ϵ, i) ∈ zip(lr_scheduler, 1:max_iter)
        # Reset gradients
        fill!(∇x_gpu, 0.0f0)
        fill!(∇y_gpu, 0.0f0)

        # Launch Metal kernel for parallel propagation & cost
        function propagation_cost_kernel(spins_gpu, cf, cost_gpu, adjs_gpu, isos_gpu)
            index = thread_index_x()
            spin = spins_gpu[index]
            isos_gpu[index] = dynamics(cf, spin)
            cost_vars = GrapeMR.cost_function(isos_gpu[index], gp.cost_function)
            adj_ini = last(cost_vars)
            adjs_gpu[index] = backward_propagation(cf, isos_gpu[index], adj_ini)
            cost_gpu[index] = first(cost_vars)
        end

        Metal.launch(propagation_cost_kernel, length(spins), (spins, cf, cost_gpu, adjs_gpu, isos_gpu))
        
        grape_output.cost_values[i, 1] = sum(cost_gpu)

        # Save final magnetization trajectory
        if i == max_iter
            append!(grape_output.isochromats, isos_gpu)
        end

        # Launch Metal kernel for gradients
        if gp.fields_opt[1] || gp.fields_opt[2]
            function gradient_kernel(adjs_gpu, dyns, ∇x_gpu, ∇y_gpu, Ix, Iy, fields_opt)
                index = thread_index_x()
                adj = adjs_gpu[index]
                dyn = dyns[index]

                if fields_opt[1]
                    atomic_add!(∇x_gpu, gradient(adj, dyn, Ix))
                end
                if fields_opt[2]
                    atomic_add!(∇y_gpu, gradient(adj, dyn, Iy))
                end
            end
            dyns_gpu = getfield.(isos_gpu, :magnetization) |> getfield.(:dynamics)
            Metal.launch(gradient_kernel, length(spins), (adjs_gpu, dyns_gpu, ∇x_gpu, ∇y_gpu, Ix, Iy, gp.fields_opt))
        end

        # Control Field update
        (u1x, u1y) = update!(cf, (∇x_gpu, ∇y_gpu), ϵ)
        cf.B1x = u1x
        cf.B1y = u1y
    end

    grape_output.control_field.B1x = u1x
    grape_output.control_field.B1y = u1y

    # Print final cost value
    final_cost = round(grape_output.cost_values[end], digits=3)
    println("\n Final Cost Function Value = $final_cost \n")

    RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


"""
    dynamics(cf::ControlField, spins::Spin)

Function that returns the Isochromat object with the already calculated dynamics.

    # Input  
    - cf::ControlField - Adjoint State
    - spins::Spin

    # Output
    - iso::Isochromat
"""
function dynamics(cf::ControlField, spin::Spin)
    mag = forward_propagation(cf, spin)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spin)
    return iso
end


"""
    gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)

# Arguments  
- χ = (::Matrix{Float64}) - Adjoint State
- M = (::Matrix{Float64}) - Forward Propagation
- H = (::Matrix) - Hamiltonian

# Outputs
- ΔJ - 1xN matrix
"""
function gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix{Int64})
    grad = zeros(Float64, 1, length(M[1,:])-1)
    for i ∈ 1:(length(M[1,:])-1)
        grad[1,i] = transpose(χ[:,i+1])*H*M[:,i+1]
    end
    return grad
end

function metal_gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix{Int64})
    grad = Metal.zeros(Float64, 1, length(M[1,:])-1; storage=Metal.SharedStorage)
    for i ∈ 1:(length(M[1,:])-1)
        grad[1,i] = transpose(χ[:,i+1])*H*M[:,i+1]
    end
    return grad
end


"""
    update(cf::ControlField, ∇xy::Tuple, ϵ::Float64)

update
# Arguments  
- cf:  (::ControlField) - Control fields struct
- ∇xy: (::Tuple) - Calculated gradients for x and y components
- ϵ:   (::Float64) - Weigth of gradient

# Outputs
- Control Field - 1xN matrix
"""
function update!(cf::ControlField, ∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ*∇xy[1]
    u1y = cf.B1y .- ϵ*∇xy[2]
    return u1x, u1y
end
