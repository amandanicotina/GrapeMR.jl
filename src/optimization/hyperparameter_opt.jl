####################################################################################
#                         Hyperparameter Optimization                              #
####################################################################################
"""
    random_sample(spins::Vector{<:Spins}, 
                  gp::GrapeParams, 
                  Tc::LinRange, 
                  max_iter::StepRange; 
                  i::Int=50, 
                  poly_start::Vector{Float64}=[1e-1, 1e-2], 
                  poly_degree::Vector{Int}=[1, 2, 3], 
                  B1ref::Float64=1.0)

Performs random sampling for hyperparameter optimization.

    ### Input
    - `spins::Vector{<:Spins}`: A vector containing spin systems for the optimization process.
    - `gp::GrapeParams`: Grape algorithm parameters, which include time points, cost function, and mask for which fields are being optimized.
    - `Tc::LinRange`: Range for the time control points for spline interpolation.
    - `max_iter::StepRange`: Range for the maximum number of iterations for the optimizer.
    - `i::Int=50`: Number of random samples to evaluate.
    - `poly_start::Vector{Float64}=[1e-1, 1e-2]`: Starting points for the polynomial learning rate.
    - `poly_degree::Vector{Int}=[1, 2, 3]`: Degrees for the polynomial learning rate.
    - `B1ref::Float64=1.0`: Reference B1 field amplitude for the RF pulse.

    ### Output
    - A random hyperparameter optimization object using random sampling.
"""
function random_sampler(spins::Vector{<:Spins}, 
            gp::GrapeParams, 
            Tc::LinRange,
            max_iter::StepRange; 
            i::Int = 50, 
            poly_start::Vector{Float64} = [1e-1, 1e-2], 
            poly_degree::Vector{Int} = [1, 2, 3], 
            B1ref::Float64 = 1.0) #, logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing))
    
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


"""
    hyperoptimization(spins::Vector{<:Spins}, 
                      gp::GrapeParams, 
                      Tc::LinRange, 
                      max_iter::Int; 
                      i::Int=5, 
                      poly_start::Vector{Float64}=[1e-1, 1e-2], 
                      poly_degree::Vector{Int}=[1, 2, 3], 
                      B1ref::Float64=1.0)

Performs hyperparameter optimization using Bayesian Optimization and Hyperband (BOHB) for selecting hyperparameters.

    ### Input
    - `spins::Vector{<:Spins}`: A vector containing spin systems for the optimization process.
    - `gp::GrapeParams`: Grape algorithm parameters, which include time points, cost function, and mask for which fields are being optimized.
    - `Tc::LinRange`: Range for the time control points for spline interpolation.
    - `max_iter::Int`: Maximum number of iterations for the optimizer.
    - `i::Int=5`: Number of optimization evaluations.
    - `poly_start::Vector{Float64}=[1e-1, 1e-2]`: Starting points for the polynomial learning rate.
    - `poly_degree::Vector{Int}=[1, 2, 3]`: Degrees for the polynomial learning rate.
    - `B1ref::Float64=1.0`: Reference B1 field amplitude for the RF pulse.

    ### Output
    - An optimized object from the BOHB optimization process, containing results and optimization history.
"""
function bohb_hyperopt(spins::Vector{<:Spins}, 
            gp::GrapeParams, 
            Tc::LinRange, 
            max_iter::Int; 
            i::Int=5, 
            poly_start::Vector{Float64} = [1e-1, 1e-2], 
            poly_degree::Vector{Int} = [1, 2, 3], 
            B1ref::Float64 = 1.0) #, logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing))
                    
    Δmax_iter   = Int(ceil(max_iter/3))
    bohb = @hyperopt for i = i, sampler = Hyperband(R=max_iter, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
        Tc = Tc,
        poly_start  = poly_start,
        poly_degree = poly_degree,
        max_iter    = range(Δmax_iter, step = Δmax_iter, stop = max_iter+Δmax_iter);
        
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





    # runs = bohb.history 
    # loss = round.(bohb.results, digits=3)  
    # c1, c2, c3, c4 = bohb.minimizer

    # configurations = [runs[i][1:end-1] for i in eachindex(runs)]
    # unique_configs = unique(configurations) 

    # """
    #     get_IDconfig_config()
    # returns a dict where the keys are the config_ids and the values
    # are the actual configurations
    # """
    # function get_IDconfig_config(unique_configs::Vector{Vector{Float64}})
    #     configID = [string("config", i) for i in eachindex(unique_configs)]
    #     dict_configs = Dict(configID[i] => unique_configs[i] for i in eachindex(unique_configs))
    #     return dict_configs
    # end


    # all_pairs = [(runs[i], loss[i]) for i in eachindex(runs)] 

    # for (keys, values) in pairs(dict_IDconfig)
    #     budgets = []
    #     unique_all_pairs = filter(x -> x[1][1:3] == dict[keys], all_pairs)
    #     budget_loss_pairs = [(unique_all_pairs[i][1][end], unique_all_pairs[i][2]) for i in 1:15]#eachindex(unique_all_pairs)]
    #     dict = Dict(keys => budget_loss_pairs)
    # end













function hband_hyperopt(spins::Vector{<:Spins}, 
            gp::GrapeParams, 
            Tc::LinRange, 
            max_iter::Int; 
            B1ref::Float64 = 1.0)

    Δmax_iter = Int(ceil(max_iter / 3))

    hb = @hyperopt for resources = max_iter, 
                       sampler = Hyperband(R = max_iter, η = 3, inner = RandomSampler()), 
                       poly_start = [1e-1, 1e-2], 
                       poly_degree = [1, 2, 3], 
                       Tc = Tc, 
                       iter = range(Δmax_iter, step = Δmax_iter, stop = max_iter + Δmax_iter)

        # Handle state if previously set
        if state !== nothing
            Tc, poly_start, poly_degree, _ = state
        end

        # Validate parameters
        if Tc >= 0.0 && poly_start >= 0.0 && poly_degree >= 1.0
            println("\n resources: ", resources, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t")

            # Generate RF control field
            control_field = spline_RF(gp.N, Tc, B1ref)

            # Set up optimization parameters
            opt_params = OptimizationParams(poly_start, poly_degree, trunc(Int, resources))
            params = Parameters(gp, opt_params)

            # Run GRAPE 
            res = grape(params, control_field, spins)

            # Extract cost from result
            cost = res.cost_values[end]
            @info "metrics" resources = resources cost = cost Tc = Tc poly_start = poly_start poly_degree = poly_degree max_iter = resources
            return cost, [Tc, poly_start, poly_degree, resources]
        else
            # Invalid parameters: return Inf for the cost, which is more informative
            return Inf, [0.0, 0.0, 0.0, 0.0]
        end
    end

    # Filter out invalid results
    hb.results = filter(x -> x != Inf, hb.results)
    hb.history = filter(x -> x != [0.0, 0.0, 0.0, 0.0], hb.history)

    return hb
end

# # Example call
# hband = @time hband_hyperopt(spins, grape_params, LinRange(0.01, 0.5, 9), 2187)



function get_resources_configurations_hband(η::Int, R::Float32)
    smax = floor(Int, log(η,R))
    B = (smax + 1)*(R)

    for s in smax:-1:0
        n = ceil(Int, (B/R)*((η^s)/(s+1)))
        r = round(R / (η^s), digits=2)
        push!(ns, n)
        print("s = $s and r = $r  and n = $n \n")
    end
end
