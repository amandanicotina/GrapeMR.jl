struct GrapeOutput
    isochromats::Vector{Isochromat}
    control_field::ControlField
    cost_values::Array{Float64}
end

function grape(op::OptimizationParams, gp::GrapeParams, cf::ControlField, spins::Vector{<:Spins})
    max_iter     = op.max_iter
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=max_iter+1) 
    cost_vals    = zeros(Float64, max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals)

    for (ϵ, i) ∈ zip(lr_scheduler, 1:max_iter)
        # ϵ   = max(ϵ, eps)
        ∇x  = zeros(Float64, 1, gp.N)
        ∇y  = zeros(Float64, 1, gp.N)
        # @show eps
        for spin ∈ spins
            # Propagation & cost
            mag = forward_propagation(cf, spin)
            dyn = GrapeMR.Magnetization(mag)
            iso = Isochromat(dyn, spin)
            grape_output.cost_values[i,1] += GrapeMR.cost_function(iso, gp.cost_function)
            cost_grad = GrapeMR.cost_function_gradient(iso, gp.cost_function)
            adj = backward_propagation(cf, iso, cost_grad)
            if i == max_iter
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


function random_sample(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::AbstractRange; i::Int = 50, poly_start::Vector{Float64} = [1e-1, 1e-2], poly_degree::Vector{Int} = [1, 2, 3])
    random_hyperopt = @hyperopt for i = i,
            Tc = Tc,
            poly_start  = poly_start,
            poly_degree = poly_degree,
            max_iter    = max_iter;
            sampler = Hyperband(R=50, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
        # BOHB
        # See: https://github.com/baggepinnen/Hyperopt.jl#bohb
        # And: https://arxiv.org/pdf/1807.01774
        # state is set by the BOHB algorithm and a KDE will estimate hyperparameters
        # that balance exploration and exploitation based on previous observations

        # print("\n", i, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t", max_iter, "\t")

        # RFs
        B1ref = 1.0
        B1x = spline_RF(gp.N, Tc)'
        B1y = spline_RF(gp.N, Tc)'
        Bz  = zeros(1,gp.N)
        control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

        # Optimize
        opt_params   = OptimizationParams(poly_start, poly_degree, max_iter)
        grape_output = grape(opt_params, gp, control_field, spins)

        cost = grape_output.cost_values[end]
        cost
    end

    return random_hyperopt
end

function hyperoptimization(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::Int; i::Int = 50, poly_start::Vector{Float64} = [1e-1, 1e-2], poly_degree::Vector{Int} = [1, 2, 3])
    bohb = @hyperopt for i = i, sampler = Hyperband(R=max_iter, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
            Tc = Tc,
            poly_start  = poly_start,
            poly_degree = poly_degree,
            max_iter = range(1, step=1, stop=max_iter);
   
        if state !== nothing
            Tc, poly_start, poly_degree, _ = state
        end
        if Tc >= 0.0 && poly_start >= 0.0 && poly_degree >= 1.0
            # print("\n resources: ", i, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t")

            # RFs
            B1ref = 1.0
            B1x = spline_RF(gp.N, Tc)'
            B1y = spline_RF(gp.N, Tc)'
            Bz  = zeros(1,gp.N)
            control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

            # Optimize
            opt_params   = OptimizationParams(poly_start, poly_degree, trunc(Int, i))
            res = grape(opt_params, gp, control_field, spins)

            cost = res.cost_values[end]
            cost, [Tc, poly_start, poly_degree, i]
        else
            1000.0, [0.0, 0.0, 0.0, 0.0]
        end
    end
    # cleanup results and history
    bohb.results = filter((x -> x != 1000.0), bohb.results)
    bohb.history = filter((x -> x != [0.0, 0.0, 0.0, 0.0]), bohb.history)

    return bohb
end