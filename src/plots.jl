function plot_magnetization(iso::Magnetization, t::Float64)
    mag  = iso.magnetization
    spin = iso.spin
    time = range(0.0, t, length = length(mag[1,:]))

    # Create a plot of the magnetization
    p = plot(time, mag[2:end,:]', label = ["Mx" "My" "Mz"], lw = 2,
        xlabel = "t [sec]",
        ylabel = "Magnitude",
        title  = "Magnetization Dynamics - Target = $(spin.target)",
        titlefontsize = 12,
        )
    return p
end

function plot_magnetization_target(iso::Magnetization)
    mag  = iso.magnetization
    spin = iso.spin
    Mxy = sqrt.(mag[2,:].^2 .+ mag[3,:].^2)
    
    p = plot(Mxy, mag[4, :], label = false, color = 1, lw = 1,
        xlims = [0.0, 1.0],
        ylims = [-1.0, 1.0],
        xlabel = "Mxy",
        ylabel = "Mz",
        title  = "Magnetization Dynamics - Target = $(spin.target)",
        titlefontsize = 12)
        scatter!([Mxy[end]], [mag[4, end]], label = false, color = 1)
    return p
end


function plot_magnetization_target3d(iso::Magnetization)
    mag  = iso.magnetization
    spin = iso.spin

    # Create a plot of the magnetization
    p = plot3d(mag[2, :], mag[3, :], mag[4, :], lw = 2,
        zlims = [-1.0, 1.0],
        xlabel = "Mx",
        ylabel = "My",
        zlabel = "Mz",
        title  = "Magnetization Dynamics - Target = $(spin.target)",
        titlefontsize = 12)
        scatter!([mag[2, end]], [mag[3, end]], [mag[4, end]])
    return p
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

function plot_control_fields(cf::ControlFields)
    if typeof(cf.B1x) == Array{Float64, 3}
        time = range(0.0, cf.t_control, length = length(cf.B1x[:,:,end]))
        Bx = cf.B1x[:,:,end]
        By = cf.B1y[:,:,end]

        p_Bx = plot(time, Bx', linewidth=2, label = false, ylabel="B1x [rads/s]", title="Optimized Control Fields", titlefontsize=12)
        p_By = plot(time, By', linewidth=2, label = false, ylabel="B1y [rads/s]", xlabel="t [s]")

        p = plot(p_Bx, p_By, layout=(2,  1))

    else

        time = range(0.0, cf.t_control, length = length(cf.B1x))
        Bx = cf.B1x
        By = cf.B1y

        p_Bx = plot(time, Bx', linewidth=2, label = false, ylabel="B1x [rads/s]", title="Initial Control Fields", titlefontsize=12)
        p_By = plot(time, By', linewidth=2, label = false, ylabel="B1y [rads/s]", xlabel="t [s]")

        p = plot(p_Bx, p_By, layout=(2,  1))
    end

    return p 
end

function plot_all_fields()
    data_x = [fields_opt.B1x[1, :, i:1000:end] for i in 1:max_iter]
    data_y = [fields_opt.B1y[1, :, i:1000:end] for i in 1:max_iter]
    p_x = plot(data_x)
    p_y = plot(data_y)

end