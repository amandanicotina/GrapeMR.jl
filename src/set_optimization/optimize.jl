struct GrapeOutput
    isochromats::Vector{Isochromat}
    control_field::ControlField
    cost_values::Array{Float64}
end

function grape(op::OptimizationParams, cf::ControlField, spins::Vector{Spin}, lr_scheduler::Poly; max_iter=2500, ϵ = 1e-4)
    cost_vals = zeros(Float64, max_iter, 1)[:]
    u1x, u1y = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals)

    for (eps, i) ∈ zip(lr_scheduler, 1:max_iter)
        ϵ   = max(ϵ, eps)
        ∇x  = zeros(Float64, 1, op.N)
        ∇y  = zeros(Float64, 1, op.N)
        @show eps
        for spin ∈ spins
            # Propagation & cost
            mag = forward_propagation(cf, spin)
            dyn = Magnetization(mag)
            iso = Isochromat(dyn, spin)
            grape_output.cost_values[i,1] += GrapeMR.cost_function(iso, op.cost_function)
            cost_grad = GrapeMR.cost_function_gradient(iso, op.cost_function)
            adj = backward_propagation(cf, iso, cost_grad)
            if i == max_iter
                push!(grape_output.isochromats, iso)
            end
            # Gradient
            if op.fields_opt[1]
                ∇x += gradient(adj, mag, Ix)
            end 
            if op.fields_opt[2]
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



