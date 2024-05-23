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

function plot_magnetization_2D(isos::Vector{Isochromat})
    p = plot(titlefontsize = 14)

    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2,:].^2 + m[3,:].^2)
        if s.target == "max" 
            plot!(p, Mxy, m[4, :], 
            label = max_label_plotted ? false : "$(s.target) - $(s.label)", color = 1, lw = 1)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = 1, markersize = 4)
            max_label_plotted = true
            
        elseif s.target == "min" 
            plot!(p, Mxy, m[4, :], 
            label = min_label_plotted ? false : "$(s.target) - $(s.label)", color = "black", lw = 1)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = "black", markersize = 4)
            min_label_plotted = true
        else
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end
    xlims!(-1.0, 1.0)
    ylims!(-1.0, 1.1)
    xlabel!("Mxy")
    ylabel!("Mz")    
    title!("Magnetization Dynamics")

    return p
end

function plot_magnetization_3D(isos::Vector{Isochromat})
    p = plot3d(titlefontsize = 14)
    
    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        if s.target == "max"
            plot3d!(p, m[2,:], m[3,:], m[4,:], 
            label = max_label_plotted ? false : "$(s.target) - $(s.label)", color = 1, lw = 1.5)
            scatter!(p, [m[2,end]], [m[3,end]], [m[4,end]], label = false, color = 1, markersize = 4)
            max_label_plotted = true
        elseif s.target == "min"
            plot3d!(p, m[2,:], m[3,:], m[4,:], 
            label = min_label_plotted ? false : "$(s.target) - $(s.label)", color = "black", lw = 1.5)
            scatter!(p, [m[2,end]], [m[3,end]], [m[4,end]], label = false, color = "black", markersize = 4)
            min_label_plotted = true
        else
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end
    xlabel!("Mx")
    ylabel!("My")
    zlabel!("Mz")
    title!("Magnetization Dynamics")

    return p
end

# Revisit this target plot functions
function plot_magnetization_target(isos::Vector{Isochromat})
    p = plot()

    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m   = iso.magnetization.dynamics
        s   = iso.spin
        Mxy = sqrt.(m[2,:].^2 .+ m[3,:].^2)
        if s.target == "max" 
            # Steady State
            ss = steady_state_matrix(s)
            Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
            plot!(p, Mxy, m[4, :], label = max_label_plotted ? false : "$(s.target) - $(s.label)", color = 1, lw = 2,
                xlims = [-1.0, 1.0],
                ylims = [-1.0, 1.1],
                xlabel = "Mxy",
                ylabel = "Mz",
                title = "Magnetization Dynamics",
                titlefontsize = 12)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = 1, markersize = 4)
            scatter!([Mxy_ss], [Mz_ss], color = 3, label = max_label_plotted ? false : "Steady State $(s.target)")
            max_label_plotted = true
            
        elseif s.target == "min" 
            # Steady State
            ss = steady_state_matrix(s)
            Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
            plot!(p, Mxy, m[4, :], label = min_label_plotted ? false : "$(s.target) - $(s.label)", color = "black", lw = 2,
                xlims = [-1.0, 1.0],
                ylims = [-1.0, 1.1],
                xlabel = "Mxy",
                ylabel = "Mz",
                title = "Magnetization Dynamics",
                titlefontsize = 12)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = "black", markersize = 4)
            scatter!([Mxy_ss], [Mz_ss], color = 4, label = min_label_plotted ? false : "Steady State $(s.target)")
            min_label_plotted = true
        else
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end

    return p
end

function plot_magnetization_target_3D(iso::Isochromat)
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
    p = plot3d(Mx, My, Mz, label = "$(s.target) - $(s.label)", color = 1, lw = 2)
        scatter!([Mx[end]], [My[end]], [Mz[end]], label = false, color = 1, markersize = 3)
    if s.target == "max"
        scatter!([Mx_ss], [My_ss], [Mz_ss], label = "Steady State", color = 3, markersize = 3)
    end
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
