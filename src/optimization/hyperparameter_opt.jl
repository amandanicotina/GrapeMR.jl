"""
    random_hyperopt(spins::Vector{<:Spins}, 
                    gp::GrapeParams, 
                    Tc::LinRange, 
                    max_iter::StepRange; 
                    i::Int = 50, 
                    poly_start::Vector{Float64} = [5e-1, 1e-1, 1e-2], 
                    poly_degree::Vector{Int} = [1, 2], 
                    B1ref::Float64 = 1.0)

Performs random sampling for hyperparameter optimization.

# Arguments
- `spins::Vector{<:Spins}`: A vector of spin systems for the optimization process.
- `gp::GrapeParams`: GRAPE algorithm parameters, including cost function and fields to optimize.
- `Tc::LinRange`: Range for time control points for spline interpolation.
- `max_iter::StepRange`: Range for maximum optimization iterations.
- `i::Int=50`: Number of random samples to evaluate.
- `poly_start::Vector{Float64}`: Initial values for polynomial learning rate.
- `poly_degree::Vector{Int}`: Degrees for the polynomial learning rate.
- `B1ref::Float64=1.0`: Reference B1 field amplitude.

# Returns
- A hyperparameter optimization object with randomly sampled configurations and their costs.
"""
function random_hyperopt(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::StepRange;
    i::Int=30,
    poly_start::Vector{Float64}=[1e-1, 2.5e-1, 5e-1, 7.5e-1],
    poly_degree::Vector{Int}=[1, 2],
    B1ref::Float64=1.0)
    #, logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing))

    random_hyperopt = @hyperopt for i = i,
        Tc = Tc,
        poly_start = poly_start,
        poly_degree = poly_degree,
        max_iter = max_iter

        # Initialize RFs
        control_field = spline_RF(gp.N, Tc, B1ref)

        # Parameters
        opt_params = OptimizationParams(poly_start, round(Int, poly_degree), trunc(Int, max_iter))
        params = Parameters(gp, opt_params)

        # Run GRAPE
        grape_output = grape(params, control_field, spins)

        # Results
        cost = grape_output.cost_values[end]
        @info "metrics" hyperopt_i = i cost = cost Tc = Tc poly_start = poly_start poly_degree = poly_degree max_iter = max_iter
        cost
    end

    return random_hyperopt

end


"""
    bohb_hyperopt(spins::Vector{<:Spins}, 
                  gp::GrapeParams, 
                  Tc::LinRange, 
                  max_iter::Int; 
                  i::Int=5, 
                  poly_start::Vector{Float64} = [5e-1, 1e-1, 1e-2], 
                  poly_degree::Vector{Int} = [1, 2], 
                  B1ref::Float64 = 1.0)

Performs hyperparameter optimization using Bayesian Optimization and Hyperband (BOHB).

# Arguments
- `spins::Vector{<:Spins}`: Vector of spin systems for the optimization process.
- `gp::GrapeParams`: GRAPE algorithm parameters, including cost function and fields to optimize.
- `Tc::LinRange`: Range for time control points for spline interpolation.
- `max_iter::Int`: Maximum number of optimization iterations.
- `i::Int=5`: Number of optimization evaluations.
- `poly_start::Vector{Float64}`: Initial values for polynomial learning rate.
- `poly_degree::Vector{Int}`: Degrees for the polynomial learning rate.
- `B1ref::Float64=1.0`: Reference B1 field amplitude.

# Returns
- A BOHB optimization object with optimized hyperparameter configurations and costs.
"""
function bohb_hyperopt(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::Int;
    i::Int=5,
    poly_start::Vector{Float64}=[5e-1, 1e-1, 1e-2],
    poly_degree::Vector{Int}=[1, 2],
    B1ref::Float64=1.0)

    # BOHB Sampler
    bohb = Hyperband(R=max_iter, η=3, 
            inner=BOHB(dims=[Hyperopt.Continuous(), 
                             Hyperopt.Continuous(), 
                             Hyperopt.Continuous(), 
                             Hyperopt.Continuous()]))

    bohb_hyperopt = @hyperopt for i = i, sampler = bohb,
        Tc = Tc,
        poly_start = poly_start,
        poly_degree = poly_degree,
        max_iter = max_iter

        # Logic for using a passed-in state
        if state !== nothing
            Tc, poly_start, poly_degree, _ = state
        end

        # Initialize RFs
        control_field = spline_RF(gp.N, Tc, B1ref)

        # Parameters
        opt_params = OptimizationParams(poly_start, round(Int, poly_degree), trunc(Int, i))
        params = Parameters(gp, opt_params)

        # Run GRAPE
        res = grape(params, control_field, spins)

        # Results
        cost = res.cost_values[end]
        @info "metrics" resources = i cost = cost Tc = Tc poly_start = poly_start poly_degree = poly_degree max_iter = i
        cost, [Tc, poly_start, poly_degree, i]
    end

    return bohb_hyperopt
end




"""
    hband_hyperopt(spins::Vector{<:Spins}, 
                   gp::GrapeParams, 
                   Tc::LinRange, 
                   max_iter::Int; 
                   B1ref::Float64 = 1.0)

Executes Hyperband optimization for selecting hyperparameters.

# Arguments
- `spins::Vector{<:Spins}`: Vector of spin systems for the optimization process.
- `gp::GrapeParams`: GRAPE algorithm parameters, including cost function and fields to optimize.
- `Tc::LinRange`: Range for time control points for spline interpolation.
- `max_iter::Int`: Maximum number of optimization iterations.
- `B1ref::Float64=1.0`: Reference B1 field amplitude.

# Returns
- A Hyperband optimization object with configurations, costs, and history.
"""
function hband_hyperopt(spins::Vector{<:Spins}, gp::GrapeParams, Tc::LinRange, max_iter::Int;
    B1ref::Float64=1.0)

    hband_sampler = Hyperband(R=max_iter, η=3, inner=RandomSampler())

    hband_hyperopt = @hyperopt for resources = max_iter, sampler = hband_sampler,
        poly_start = [1e-1, 1e-2],
        poly_degree = [1, 2],
        Tc = Tc

        # Handle state if previously set
        if state !== nothing
            Tc, poly_start, poly_degree = state
        end

        # Initialize RF 
        control_field = spline_RF(gp.N, Tc, B1ref)

        # Parameters
        opt_params = OptimizationParams(poly_start, round(Int, poly_degree), trunc(Int, resources))
        params = Parameters(gp, opt_params)

        # Run GRAPE 
        res = grape(params, control_field, spins)

        # Results
        cost = res.cost_values[end]
        @info "metrics" resources = resources cost = cost Tc = Tc poly_start = poly_start poly_degree = poly_degree max_iter = resources
        return cost, [Tc, poly_start, poly_degree]
    end

    return hband_hyperopt
end






"""
    get_resources_configurations_hband(η::Int, R::Float64)

Calculates the resources and configurations for each stage of the Hyperband algorithm.

# Arguments
- `η::Int`: Downsampling factor, typically used in Hyperband.
- `R::Float64`: Maximum resource budget (e.g., number of iterations).

# Returns
- Prints the configurations `n` and resources `r` for each Hyperband outer loop.
"""
function get_resources_configurations_hband(η::Int, R::Real)
    smax = floor(Int, log(η, R))
    B = (smax + 1) * (R)
    ns = []
    for s in smax:-1:0
        n = ceil(Int, (B / R) * ((η^s) / (s + 1)))
        r = round(R / (η^s), digits=2)
        push!(ns, n)
        print("s = $s and r = $r  and n = $n \n")
    end
end

# R = 1500
# get_resources_configurations_hband(3, R)

# Spearman's Rank

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
