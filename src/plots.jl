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

function plot_magnetization_target(iso::Isochromat)
    mag  = iso.magnetization.dynamics
    spin = iso.spin
    Mxy = sqrt.(mag[2,:].^2 .+ mag[3,:].^2)

    Mz_tar = [0.0]
    Mxy_tar = sqrt((0.5)^2 + (0.5)^2)
    p = plot(Mxy, mag[4, :], label = false, color = 1, lw = 1,
        xlims = [-0.001, 1.0],
        ylims = [-1.0, 1.0],
        xlabel = "Mxy",
        ylabel = "Mz",
        title  = "Magnetization Dynamics - $(spin.label)",
        titlefontsize = 12)
        scatter!([Mxy[end]], [mag[4, end]], label = false, color = 1)
        scatter!([Mxy_tar], [Mz_tar], label = false, color = 2)
    return p
end


function plot_magnetization_Mz_Mxy(isos::Vector{Isochromat})
    p = plot()
    #curve_color = range(1.0, length(grape_output.isochromats))

    for (idx, iso) in enumerate(isos)
        mag = iso.magnetization.dynamics
        spin = iso.spin
        Mxy = sqrt.(mag[2,:].^2 .+ mag[3,:].^2)

        plot!(p, Mxy, mag[4, :], label = "$(spin.target)", color = idx, lw = 2,
            xlims = [-1.0, 1.0],
            ylims = [-1.0, 1.0],
            xlabel = "Mxy",
            ylabel = "Mz",
            title = "Magnetization Dynamics - $(spin.label)",
            titlefontsize = 12)
        
        scatter!(p, [Mxy[end]], [mag[4, end]], label = false, color = idx, markersize = 4)
    end

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

function plot_cost_values(cost::Vector{Float64}, op::OptimizationParams)
    p = plot(cost, label = op.cost_function, lw = 2,
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

    p_Bx = plot(time, Bx', linewidth=2, label = false, ylabel="u1x", title="Control Fields", titlefontsize=12)
    p_By = plot(time, By', linewidth=2, label = false, ylabel="u1y", xlabel="t [s]")

    p = plot(p_Bx, p_By, layout=(2,  1))

    return p 
end
