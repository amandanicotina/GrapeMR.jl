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
    max_iter     = op.max_iter
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=max_iter+1) 
    cost_vals    = zeros(Float64, max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
    ∇x  = zeros(Float64, 1, gp.N)
    ∇y  = zeros(Float64, 1, gp.N)
    
    for (ϵ, i) ∈ zip(lr_scheduler, 1:max_iter)
        # ϵ   = max(ϵ, eps)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)

        # Propagation & cost
        iso = dynamics.(cf, spins)
        cost_vars = GrapeMR.cost_function.(iso, gp.cost_function)
        cost_val = sum(first.(cost_vars)) 
        adj_ini  = last.(cost_vars)    
        adj      = backward_propagation.(Ref(cf), iso, adj_ini)
        grape_output.cost_values[i,1] = cost_val

        # Save final magnetization trajectory
        if i == max_iter
            append!(grape_output.isochromats, iso)
        end
        # Gradient
        if gp.fields_opt[1]
            ∇x = sum(gradient.(adj, getfield.(getfield.(iso, :magnetization), :dynamics), Ref(Ix)))
        end 
        if gp.fields_opt[2]
            ∇y = sum(gradient.(adj, getfield.(getfield.(iso, :magnetization), :dynamics), Ref(Iy)))
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
    dynamics(cf::ControlField, spins::Vector{Spin})

Function that returns the Isochromat object with the already calculated dynamics.

    # Input  
    - cf::ControlField - Adjoint State
    - spins::Vector{Spin}

    # Output
    - iso::Vector{Isochromat}
"""
function dynamics(cf::ControlField, spins::Spin)
    mag = forward_propagation(cf, spins)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spins)
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







function random_sample(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::StepRange; i::Int = 50, poly_start::Vector{Float64} = [1e-1, 1e-2], poly_degree::Vector{Int} = [1, 2, 3], B1ref::Float64 = 1.0)#, logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing))
    random_hyperopt = @hyperopt for i = i,
            Tc = Tc,
            poly_start  = poly_start,
            poly_degree = poly_degree,
            max_iter    = max_iter
        # RFs
        control_field = spline_RF(gp.N, Tc, B1ref);

        # Optimize
        opt_params   = OptimizationParams(poly_start, poly_degree, max_iter);
        params       = Parameters(gp, opt_params)
        grape_output = grape(params, control_field, spins);

        cost = grape_output.cost_values[end];
        @info "metrics" hyperopt_i=i cost=cost Tc=Tc poly_start=poly_start poly_degree=poly_degree max_iter=max_iter
        cost
    end

    return random_hyperopt
end

function hyperoptimization(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::Int; i::Int=5, poly_start::Vector{Float64} = [1e-1, 1e-2], poly_degree::Vector{Int} = [1, 2, 3], B1ref::Float64 = 1.0)#, logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing))

    bohb = @hyperopt for i = i, sampler = Hyperband(R=max_iter, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
            Tc = Tc,
            poly_start  = poly_start,
            poly_degree = poly_degree,
            max_iter = range(1, step=1, stop=max_iter);

        # we're cheating here a bit, we're not using the sampled max_iter
        # this is a workaround to have the max_iter/budget in the history
        if state !== nothing
            Tc, poly_start, poly_degree, _ = state
        end
        if Tc >= 0.0 && poly_start >= 0.0 && poly_degree >= 1.0
             print("\n resources: ", i, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t")

            # RFs
            control_field = spline_RF(gp.N, Tc, B1ref);

            # Optimize
            opt_params = OptimizationParams(poly_start, poly_degree, trunc(Int, i))
            params     = Parameters(gp, opt_params)
            res = grape(params, control_field, spins)

            cost = res.cost_values[end]
            @info "metrics" resources=i cost=cost Tc=Tc poly_start=poly_start poly_degree=poly_degree max_iter=i
            cost, [Tc, poly_start, poly_degree, i]
        else
            1000.0, [0.0, 0.0, 0.0, 0.0]
        end
    end
    # cleanup results and history
    bohb.results = filter((x -> x != 1000.0), bohb.results)
    bohb.history = filter((x -> x != [0.0, 0.0, 0.0, 0.0]), bohb.history)

    # close(logger)
    return bohb
end

