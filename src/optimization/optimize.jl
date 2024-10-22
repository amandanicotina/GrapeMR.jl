const wandb_project::String = "GrapeMR"
Base.broadcastable(cf::ControlField) = Ref(cf)


struct GrapeOutput
    isochromats::Vector{Isochromat}
    control_field::ControlField
    cost_values::Array{Float64}
    params::Parameters
end


"""
    grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})

Implements the grape algorithm [khaneja2005optimal](@cite) 

    ### Input
    - `op::OptimizationParams`: Parameters for the optimization itself: max iterations, 
    - `gp::GrapeParams`: Parameters related to Grape itself: time points, cost function, mask for which fields are being optimized.
    - `cf::ControlField`: Initial control field - spline function -
    - `spins::Vector{<:Spins}`: Vector with all spins included in the optimization

    ### Output
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

gradient
    # Input  
    - χ = (::Matrix{Float64}) - Adjoint State
    - M = (::Matrix{Float64}) - Forward Propagation
    - H = (::Matrix) - Hamiltonian

    # Output
    - ΔJ - 1xN matrix
"""
function gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix{Int64})
    grad = zeros(Float64, 1, length(M[1,:])-1)
    for i ∈ 1:(length(M[1,:])-1)
        grad[1,i] = transpose(χ[:,i+1])*H*M[:,i+1]
    end
    return grad
end


"""
    update(cf::ControlField, ∇xy::Tuple, ϵ::Float64)

update
    # Input  
    - cf:  (::ControlField) - Control fields struct
    - ∇xy: (::Tuple) - Calculated gradients for x and y components
    - ϵ:   (::Float64) - Weigth of gradient

    # Output
    - Control Field - 1xN matrix
"""
function update!(cf::ControlField, ∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ*∇xy[1]
    u1y = cf.B1y .- ϵ*∇xy[2]
    return u1x, u1y
end
