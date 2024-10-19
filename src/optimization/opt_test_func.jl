## Original grape function
function old_grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    max_iter     = op.max_iter
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=max_iter+1) 
    cost_vals    = zeros(Float64, max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
    ∇xoutput  = Vector{Matrix{Float64}}()
    ∇youtput  = Vector{Matrix{Float64}}()
    
    for (ϵ, i) ∈ zip(lr_scheduler, 1:max_iter)
        # ϵ   = max(ϵ, eps)
        ∇x  = zeros(Float64, 1, gp.N)
        ∇y  = zeros(Float64, 1, gp.N)
        for spin ∈ spins
            # Propagation & cost
            mag = forward_propagation(cf, spin)
            dyn = GrapeMR.Magnetization(mag)
            iso = Isochromat(dyn, spin)
            (cost, cost_grad) = GrapeMR.cost_function(iso, gp.cost_function)
            grape_output.cost_values[i,1] += cost
            # cost_grad = GrapeMR.cost_function_gradient(iso, gp.cost_function)
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
            push!(∇xoutput, ∇x)
            push!(∇youtput, ∇y)
        end

        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x = u1x
        cf.B1y = u1y
    end

    grape_output.control_field.B1x = u1x
    grape_output.control_field.B1y = u1y

    RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


function no_threads_grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    max_iter     = op.max_iter
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=max_iter+1) 
    cost_vals    = zeros(Float64, max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
   
    ∇x   = zeros(Float64, 1, gp.N)
    ∇y   = zeros(Float64, 1, gp.N)
    isos = Vector{Isochromat}(undef, length(spins))
    adjs = Vector{Matrix{Float64}}(undef, length(spins))

    for (ϵ, i) ∈ zip(lr_scheduler, 1:max_iter)
        # ϵ   = max(ϵ, eps)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)

        # Propagation & cost
        for (j, spin) ∈ enumerate(spins)
            isos[j]   = dynamics(cf, spin)
            cost_vars = GrapeMR.cost_function(isos[j], gp.cost_function)
            cost_val  = first(cost_vars)
            adj_ini   = last(cost_vars)    
            adjs[j]   = backward_propagation(cf, isos[j], adj_ini)          
            grape_output.cost_values[i, 1] += cost_val
        end
        # Save final magnetization trajectory
        if i == max_iter
            append!(grape_output.isochromats, isos)
        end
        # Gradient
        if gp.fields_opt[1]
            ∇x = sum(gradient.(adjs, getfield.(getfield.(isos, :magnetization), :dynamics), Ref(Ix)))
        end 
        if gp.fields_opt[2]
            ∇y = sum(gradient.(adjs, getfield.(getfield.(isos, :magnetization), :dynamics), Ref(Iy)))
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
    

function threads_grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    max_iter     = op.max_iter
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=max_iter+1) 
    cost_vals    = zeros(Float64, max_iter, 1)[:]
    u1x, u1y     = [], []
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
   
    ∇x   = zeros(Float64, 1, gp.N)
    ∇y   = zeros(Float64, 1, gp.N)
    isos = Vector{Isochromat}(undef, length(spins))
    adjs = Vector{Matrix{Float64}}(undef, length(spins))

    for (ϵ, i) ∈ zip(lr_scheduler, 1:max_iter)
        # ϵ   = max(ϵ, eps)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)

        # Propagation & cost
        tmp_costs = zeros(Float64, length(spins))
        Threads.@threads for j in 1:length(spins)
            spin = spins[j]
            isos[j]      = dynamics(cf, spin)
            cost_vars    = GrapeMR.cost_function(isos[j], gp.cost_function)
            adj_ini      = last(cost_vars)    
            adjs[j]      = backward_propagation(cf, isos[j], adj_ini)          
            tmp_costs[j] = first(cost_vars)
        end
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
            ∇y = sum(gradient.(adjs, dyns, Ref(Iy)))
        end 
        
        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x = u1x
        cf.B1y = u1y
    end

    grape_output.control_field.B1x = u1x
    grape_output.control_field.B1y = u1y
    
    # Print final cost value
    final_cost = round(grape_output.cost_values[end], digits = 3)
    println("\n Final Cost Function Value = $final_cost \n")
    
    RF_pulse_analysis(grape_output.control_field)

    return grape_output
end
    