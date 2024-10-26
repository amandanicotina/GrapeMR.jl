####################################################################################
#                         Hyperparameter Optimization                              #
####################################################################################
"""
    random_sampler(spins::Vector{<:Spins}, 
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
    bohb_hyperopt(spins::Vector{<:Spins}, 
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





# function get_budget_loss_pair(hpo)
#     runs = hpo.history
#     loss = hpo.results    
#     c1, c2, c3, c4 = hpo.minimizer

#     # Gets random number from a range on integers
#     get_random_configuration(min::Int, max::Int) = rand(min:max)

#     # Creates tuples
#     budget_loss_pairs = [(round(runs[i][end], digits=2), loss[i]) for i in eachindex(runs)]
#     config_loss_pairs = [(runs[i], loss[i]) for i in eachindex(runs)]
#     configurations = [runs[i][1:end-1] for i in eachindex(runs)]
    
#     # Get unique arrays
#     unique_budgets = unique([pair[1] for pair in budget_loss_pairs])[end-2:end]
#     unique_configs = unique(configurations) # Target
#     unique_configs_budgets = [filter(x -> x[1][1:3] == unique_config, config_loss_pairs) for unique_config in unique_configs] 

   
#     # Chose configuration that minimizes the optimization
#     minimizer_configs = filter(x -> isapprox(x[1][1:3], [c1, c2, c3], atol=1e-8), vcat(unique_configs_budgets...))
#     loss_values = [x[2] for x in minimizer_configs]



#     # Dict for each unique configuration and a budget
#     unique_configs = Dict()
#     for (config, value) in minimizer_configs
#         key = config[end]  # Use the last element of the config as the key
#         if !haskey(unique_configs, key) || unique_configs[key][2] < value
#             unique_configs[key] = (config, value)  # Store the best value for each key
#         end
#     end



#     # Creates Dict with all losses related to each budget
#     dict_budget_loss = Dict(budget => [pair[2] for pair in budget_loss_pairs if pair[1] == budget] for budget in unique_budgets)

#     return dict_budget_loss
# end

# filtered_dict = sort(filter(kv -> kv[1] > 200, unique_configs))

# spearman_matrix = fill(undef, length(filtered_dict), length(filtered_dict))  # Initialize a matrix with NaN values


# # Calculate Spearman correlation for each combination of budgets
# for i in 1:n_budgets
#     for j in i:n_budgets
#         if i != j
#             values_i = budget_groups[budgets[i]]
#             values_j = budget_groups[budgets[j]]
#             n_i = length(values_i)
#             n_j = length(values_j)

#             if n_i == n_j && n_i > 1  # Ensure the sizes match
#                 ranked_i = rank(values_i)
#                 ranked_j = rank(values_j)
#                 d_squared_sum = sum((ranked_i[k] - ranked_j[k])^2 for k in 1:n_i)
#                 spearman_corr = 1 - (6 * d_squared_sum) / (n_i * (n_i^2 - 1))
#                 spearman_matrix[i, j] = spearman_corr
#                 spearman_matrix[j, i] = spearman_corr  # The matrix is symmetric
#             end
#         end
#     end
# end

# # Plot the heatmap
# heatmap(spearman_matrix, title="Spearman Rank Correlation Across Budgets", xlabel="Budgets", ylabel="Budgets", color=:inferno)


# using Hyperopt

# function hband_hyperopt(spins::Vector{<:Spins}, 
#             gp::GrapeParams, 
#             Tc::LinRange, 
#             max_iter::Int; 
#             B1ref::Float64 = 1.0)

#     Δmax_iter = Int(ceil(max_iter / 3))

#     hb = @hyperopt for resources = max_iter, 
#                        sampler = Hyperband(R = max_iter, η = 3, inner = RandomSampler()), 
#                        poly_start = [1e-1, 1e-2], 
#                        poly_degree = [1, 2, 3], 
#                        Tc = Tc, 
#                        iter = range(Δmax_iter, step = Δmax_iter, stop = max_iter + Δmax_iter)

#         # Handle state if previously set
#         if state !== nothing
#             Tc, poly_start, poly_degree, _ = state
#         end

#         # Validate parameters
#         if Tc >= 0.0 && poly_start >= 0.0 && poly_degree >= 1.0
#             println("\n resources: ", resources, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t")

#             # Generate RF control field
#             control_field = spline_RF(gp.N, Tc, B1ref)

#             # Set up optimization parameters
#             opt_params = OptimizationParams(poly_start, poly_degree, trunc(Int, resources))
#             params = Parameters(gp, opt_params)

#             # Run GRAPE 
#             res = grape(params, control_field, spins)

#             # Extract cost from result
#             cost = res.cost_values[end]
#             @info "metrics" resources = resources cost = cost Tc = Tc poly_start = poly_start poly_degree = poly_degree max_iter = resources
#             return cost, [Tc, poly_start, poly_degree, resources]
#         else
#             # Invalid parameters: return Inf for the cost, which is more informative
#             return Inf, [0.0, 0.0, 0.0, 0.0]
#         end
#     end

#     # Filter out invalid results
#     hb.results = filter(x -> x != Inf, hb.results)
#     hb.history = filter(x -> x != [0.0, 0.0, 0.0, 0.0], hb.history)

#     return hb
# end

# # Example call
# hband = @time hband_hyperopt(spins, grape_params, LinRange(0.01, 0.5, 9), 2187)



η = 3
R = 3^7
smax = floor(Int, log(η,R))

B = (smax + 1)*(R)

ns=[]
for s in smax:-1:0
    n = ceil(Int, (B/R)*((η^s)/(s+1)))
    r = round(R / (η^s), digits=2)
    push!(ns, n)
    print("s = $s and r = $r  and n = $n \n")
end

sn = sum(ns)