function plot_magnetization(iso::Magnetization, t::Float64)
    mag  = iso.magnetization[1]
    spin = iso.spin[1]
    time = range(0.0, t, length = length(mag[1,:]))

    # Create a plot of the magnetization
    p = plot(time, mag[2:end,:]', label = ["Mx" "My" "Mz"], lw = 2,
        xlabel = "t [sec]",
        ylabel = "Magnitude",
        title  = "Magnetization Dynamics - Target = $(spin.target)",
        titlefontsize = 12,
        )
    

    # Return the array of plots
    return p
end

function plot_magnetization_target()
    
end

function plot_cost_values(cost::Array)
    p = plot(cost', label = "Euclidean Distance", lw = 2,
    xlabel = "Iterations",
    ylabel = "Cost Value",
    title  = "Cost Function Convergence",
    titlefontsize = 12,
    )

    return p
end

function plot_control_fields()
    
end