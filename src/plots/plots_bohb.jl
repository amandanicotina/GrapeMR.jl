


function plot_bohb(bohb)
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