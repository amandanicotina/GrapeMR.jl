# TODO multiple dispatch if input is ::Vector{Isochromat} or ::GrapeOutput

function plot_magnetization_time(iso::Isochromat, t::Float64)
    mag  = iso.magnetization.dynamics
    spin = iso.spin
    time = range(0.0, t, length = length(mag[1,:]))

    # Create a plot of the magnetization
    p = plot(time, mag[2:end,:]', label = ["Mx" "My" "Mz"], lw = 2,
        xlabel = "t [sec]",
        ylabel = "Magnitude",
        title  = "Magnetization Dynamics - Sample = $(spin.label), Target = $(spin.target)",
        titlefontsize = 12,
        )
    return p
end



function plot_magnetization_Mz_Mxy(isos::Vector{Isochromat})
    p = plot()

    max_label_plotted = false
    min_label_plotted = false

    for iso in isos
        mag = iso.magnetization.dynamics
        spin = iso.spin
        Mxy = sqrt.(mag[2,:].^2 .+ mag[3,:].^2)
        if spin.target == "max" 
            plot!(p, Mxy, mag[4, :], label = max_label_plotted ? false : "$(spin.target) - $(spin.label)", color = 1, lw = 1.5,
                xlims = [-1.0, 1.0],
                ylims = [-1.0, 1.1],
                xlabel = "Mxy",
                ylabel = "Mz",
                title = "Magnetization Dynamics",
                titlefontsize = 12)
            scatter!(p, [Mxy[end]], [mag[4, end]], label = false, color = 1, markersize = 4)
            max_label_plotted = true
            
        elseif spin.target == "min" 
            plot!(p, Mxy, mag[4, :], label = min_label_plotted ? false : "$(spin.target) - $(spin.label)", color = "black", lw = 1.5,
                xlims = [-1.0, 1.0],
                ylims = [-1.0, 1.1],
                xlabel = "Mxy",
                ylabel = "Mz",
                title = "Magnetization Dynamics",
                titlefontsize = 12)
            scatter!(p, [Mxy[end]], [mag[4, end]], label = false, color = "black", markersize = 4)
            min_label_plotted = true
        end
 
    end

    return p
end

function plot_magnetization_target(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics

    # Steady State
    ss = steady_state_matrix(s)
    Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z

    # Magnetization
    Mxy = sqrt.((m[2,:]).^2 + (m[3,:]).^2)
    Mz  = m[4,:]

    p = plot(Mxy, Mz, label = false, color = 1, lw = 1,
        xlims = [-0.001, 1.0],
        ylims = [-1.0, 1.0],
        xlabel = "Mxy",
        ylabel = "Mz",
        title  = "Magnetization Dynamics - $(s.label)",
        titlefontsize = 12)
        scatter!([Mxy[end]], [Mz[end]], label = false, color = 1)
        scatter!([Mxy_ss], [Mz_ss], color = 3, label = "Steady State")
    return p
end

function plot_magnetization_target_3d(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics

    # Steady State
    ss = steady_state_matrix(s)
    Mx_ss, My_ss, Mz_ss = ss.x, ss.y, ss.z

    # Magnetization
    Mx = m[2, :]
    My = m[3, :]
    Mz = m[4, :]

    # Plot
    p = plot3d(Mx, My, Mz, label = "Dynamics", color = 1)
    scatter!([Mx[end]], [My[end]], [Mz[end]], label = false, color = 1, markersize = 3)
    scatter!([Mx_ss], [My_ss], [Mz_ss], label = "Steady State", color = 3, markersize = 3)

    xlabel!("Mx")
    ylabel!("My")
    zlabel!("Mz")
    title!("Magnetization Dynamics")

    return p
end



function plot_cost_values(cost::Vector{Float64}, op::OptimizationParams)
    p = plot(cost, label = string(op.cost_function), lw = 2,
    xlabel = "Iterations",
    ylabel = "Cost Value",
    title  = "Cost Function Convergence",
    titlefontsize = 12,
    )

    return p
end

function plot_control_fields(cf::ControlField)
    time = range(0.0, cf.t_control, length = length(cf.B1x))
    Bx = cf.B1x
    By = cf.B1y

    p_Bx = plot(time, Bx', linewidth=2, label = false, ylabel="B1x [Hz]", title="Control Fields", titlefontsize=12)
    p_By = plot(time, By', linewidth=2, label = false, ylabel="B1y [Hz]", xlabel="t [s]")

    p = plot(p_Bx, p_By, layout=(2,  1))

    return p 
end
