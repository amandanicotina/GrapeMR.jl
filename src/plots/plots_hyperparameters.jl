function get_hyperopt_parameters(ho::Hyperoptimizer)
    ct = [history[1] for history in ho.history]
    ps = [history[2] for history in ho.history]
    pd = [history[3] for history in ho.history]
    it = [history[4] for history in ho.history]
    return ct, ps, pd, it
end


function plot_hyperopt_history(ho::Hyperoptimizer; title::String = "")
    trials = 1:Int(length(ho.results))
    cost = ho.results
    cost_min = round(ho.minimum, digits=3)
    colors = color_palette(10)
    
    pHistory = initialize_plot("Optimization History - $title", "Trials", "Cost Value")
    
    plot!(pHistory, trials, cost, label = false, lw = 3.5, color = colors[2], alpha = 0.7)
    scatter!(pHistory, trials, cost, label = false, color = colors[2], marker = :dot, markersize = 5, alpha = 1.0)
    
    #ylims!(pHistory, -0.01, 1.0)
    annotate!(pHistory, ((length(ho.results)/2), (1-2*cost_min), text("Cost min = $cost_min", :black, 15)))

    return pHistory
end


function plot_cost_grape_runs(ho::Hyperoptimizer; plotlog::Bool=false, title::String = "")
    trials = collect(1:Int(length(ho.results)))
    t = Int.(round.(trials./100))
    # Hyperparameters
    _, _, _, it = get_hyperopt_parameters(ho)
    cost = ho.results
    cost_min = round(ho.minimum, digits = 3)
    it_min = ho.minimizer[4]
    if plotlog
        pHistory = initialize_plot(title, "log(Grape Runs)", "Cost Value")
        scatter!(pHistory, log.(it), cost, label = false, markersize = 6, color = :viridis, zcolor=t, colorbar_title = "Trials 10²") #")#
        scatter!(pHistory, [log(it_min)], [cost_min], label = "Min = $cost_min", 
                    marker = :star5, markersize = 10, color = :red)
        #ylims!(-0.02, 1.0)
    else
        pHistory = initialize_plot(title, "Grape Runs", "Cost Value")
        scatter!(pHistory, it, cost, label = false, markersize = 6, color = :viridis, zcolor=t,colorbar_title = "Trials 10²")
        scatter!(pHistory, [it_min], [cost_min], label = "Min = $cost_min", 
                    marker = :star5, markersize = 10, color = :red)
        #|ylims!(-0.02, 1.0)
    end
    return pHistory
end

function plot_cost_hyperparam(ho::Hyperoptimizer; grid_plot::Bool=true, title::String="")
    ct, ps, pd, it = get_hyperopt_parameters(ho)

    trials = 1:Int(length(ho.results))
    cost = ho.results
    cost_min = round(ho.minimum, digits = 3)
    ct_min = ho.minimizer[1]
    ps_min = ho.minimizer[2]
    pd_min = ho.minimizer[3]
    it_min = ho.minimizer[4]

    if grid_plot
        p1 = initialize_plot("", "Control Time", "Cost Value")
        p2 = initialize_plot("", "Poly Start", "Cost Value")
        p3 = initialize_plot("", "Poly Degree", "Cost Value")
        p4 = initialize_plot("", "Grape Runs (10²)", "Cost Value")
    
        scatter!(p1, ct, cost, color=:viridis, label = false, zcolor=trials)
        scatter!(p1, [ct_min], [cost_min], label = false, 
            marker = :star5, markersize = 7, color = :red)
        #ylims!(-0.02, 1.0)

        scatter!(p2, round.(ps, digits=1), cost, color=:viridis, label = false, zcolor=trials, colorbar_title = "Trials")
        scatter!(p2, [ps_min], [cost_min], label = false, 
            marker = :star5, markersize = 7, color = :red)
        #ylims!(-0.02, 1.0)

        scatter!(p3, pd, cost, color=:viridis, label = false, zcolor=trials)
        scatter!(p3, [pd_min], [cost_min], label = false, 
            marker = :star5, markersize = 7, color = :red)
        #ylims!(-0.02, 1.0)

        scatter!(p4, it*1e-2, cost, color=:viridis, label = false, zcolor=trials, colorbar_title = "Trials")
        scatter!(p4, [it_min*1e-2], [cost_min], label = false, 
            marker = :star5, markersize = 7, color = :red)
        #ylims!(-0.02, 1.0)

        pScatter1 = plot(p1, p2, p3, p4, layout = (2,2))
        return pScatter1

    else
        p1 = initialize_plot("", "Control Time", "Cost Value")
        p2 = initialize_plot("", "Poly Start", "")
        p3 = initialize_plot("", "Poly Degree", "Cost Value")
        p4 = initialize_plot("", "Grape Runs (10²)", "")
    
        scatter!(p1, ct, cost, color=:viridis, label = false, zcolor=trials, colorbar=false)
        scatter!(p1, [ct_min], [cost_min], label = false, marker = :star5, markersize = 8, color = :red)
        annotate!(p1, (ct_min, (1-cost_min), text("Cost min = $cost_min", :black, 15)))
        #ylims!(-0.02, 1.0)

        scatter!(p2, round.(ps, digits=1), cost, color=:viridis, label = false, zcolor=trials)
        scatter!(p2, [ps_min], [cost_min], label = false, marker = :star5, markersize = 8, color = :red)
        #ylims!(-0.02, 1.0)

        scatter!(p3, pd, cost, color=:viridis, label = false, zcolor=trials, colorbar=false)
        scatter!(p3, [pd_min], [cost_min], label = false, marker = :star5, markersize = 8, color = :red)
        annotate!(p3, ((1.5), (1-cost_min), text("Cost min = $cost_min", :black, 15)))
        #ylims!(-0.02, 1.0)

        scatter!(p4, it*1e-2, cost, color=:viridis, label = false, zcolor=trials)
        scatter!(p4, [it_min*1e-2], [cost_min], label = false, 
            marker = :star5, markersize = 8, color = :red)
        #ylims!(-0.02, 1.0)

        pScatter2 = plot(p1, p2, title = title, layout = (1,2))
        pScatter3 = plot(p3, p4, title = title, layout = (1,2))
        display(pScatter2), display(pScatter3)
        return pScatter2, pScatter3
    end
end



function plot_hyperopt_contour(ho::Hyperoptimizer)
    # Result Cost
    cost = [results[1] for results in ho.results]

    # Hyperparameters
    ct, ps, pd, it = get_hyperopt_parameters(ho)

    p1 = initialize_plot("", "Grape Runs", "Control Time")
    p2 = initialize_plot("", "Poly Degree", "Poly Start")
    p3 = initialize_plot("", "Poly Start", "Control Time")
    p4 = initialize_plot("", "Grape Runs", "Poly Start")

    contourf!(p1, cost, color=:viridis)
    contourf!(p2, pd, ps, cost, color=:viridis)
    contourf!(p3, ps, ct, cost, color=:viridis)
    contourf!(p4, it, ps, cost, color=:viridis)

    pCountour = plot(p1, p2, p3, p4, layout = (2,2))
    
    return pCountour
end

function plot_bohb(bohb::Hyperoptimizer)
    cost  = []
    t_c   = []
    start = []
    deg   = []
    iter   = []
    for i ∈ eachindex(bohb.results)
        push!(cost, bohb.results[i][1])
        push!(t_c, bohb.history[i][1])
        push!(start, bohb.history[i][2])
        push!(deg, bohb.history[i][3])
        push!(iter, bohb.history[i][4])
    end

    time_opt = bohb.minimizer[1]
    iter_opt = bohb.minimizer[4]
    cost_opt = round(bohb.minimum, digits=4)
    order = collect(1:length(t_c))

    # Time
    p_cost = scatter(t_c, iter, zcolor=cost, 
        markerstrokecolor = :auto, label = false,
        xlabel = "Control time [s]", ylabel = "Iterations ", colorbar_title = "Cost Value",
        title = "Hyperparameter Tuning - Control Time",
        color = :viridis)
        scatter!([time_opt], [iter_opt], label = "Minimum = $cost_opt", 
        marker = :star5, markersize = 8, color = :red)

    # Iter
    p_order =  scatter(cost, log.(iter), zcolor=order, 
        markerstrokecolor = :auto, label = false, 
        xlabel = "Cost Function", ylabel = "Iterations - log scale", colorbar_title = "Order",
        title = "Hyperparameter Tuning - Sampling Order", 
        color = :viridis)
        scatter!([cost_opt], [log.(iter_opt)], label = "Minimum = $cost_opt", 
        marker = :star5, markersize = 8, color = :red)
        
        display(p_cost)
        display(p_order)

    return p_cost, p_order
end

# p1, p2 = bohb_params(rand_hopt)



function plot_evaluations(bohb)
    plots = []
    for (i, i_param) in enumerate(bohb.params)
        for (j, j_param) in enumerate(bohb.params)
            if i == j
                hist = histogram([val[i] for val in bohb.history], bins=10, xlabel=i_param, ylabel="Number of samples", legend=false)
                push!(plots, hist)
            elseif j > i
                # we don't care about these
                continue
            else
                # plot x,y with colour equivalent to cost value
                xs = [val[j] for val in bohb.history]
                ys = [val[i] for val in bohb.history]
                xy = scatter(
                    xs, ys, zcolor=bohb.results, xlabel=j_param, ylabel=i_param, legend=false,
                    color = :viridis
                )
                scatter!([bohb.minimizer[j]], [bohb.minimizer[i]], marker = :star5, markersize = 8, color = :red)
                push!(plots, xy)
            end
        end
    end
    grid_plot = plot(plots..., layout=@layout([° _ _ _; ° ° _ _; ° ° ° _; ° ° ° °]))
    
    return grid_plot
end