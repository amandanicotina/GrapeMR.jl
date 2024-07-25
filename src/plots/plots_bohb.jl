
function plot_bohb(bohb)
    cost  = []
    t_c   = []
    start = []
    deg   = []
    max   = []
    for i ∈ eachindex(bohb.results)
        push!(cost, bohb.results[i][1])
        push!(t_c, bohb.history[i][1])
        push!(start, bohb.history[i][2])
        push!(deg, bohb.history[i][3])
        push!(max, bohb.history[i][4])
    end
    function get_eps(start, deg, max_iter)
        return start / (1 - (max_iter/2 - 1) / max_iter)^deg
    end
    eps = [get_eps(start, deg, it) for (start, deg, it) in zip(start, deg, max)]
    t_min = bohb.minimizer[1]
    st_min = bohb.minimizer[2]
    d_min = bohb.minimizer[3]
    m_min = bohb.minimizer[4]
    c_min = round(bohb.minimum[1], digits=4)
    order = collect(1:length(t_c))
    m_iter = 1000:100:5000
    p_cost =  scatter(t_c, max, zcolor=cost, 
        markerstrokecolor = :auto, label = false,
        xlabel = "Control time [s]", ylabel = "Iter", colorbar_title="Cost Value",
        title = "Optimization",
        color = :viridis)
        scatter!([t_min], [m_min], label = "Minimum = $c_min", 
        marker = :star5, markersize = 8, color = :red)
    p_order =  scatter(cost, m_iter, zcolor=order, 
        markerstrokecolor = :auto, label = false, 
        xlabel = "Cost Function", ylabel = "Resources - log scale", colorbar_title="Iterations",
        title = "Hyperparameter Tuning - Random Sampler", 
        color = :viridis)
        scatter!([c_min], [m_min], label = "Minimum = $c_min", 
        marker = :star5, markersize = 8, color = :red)

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