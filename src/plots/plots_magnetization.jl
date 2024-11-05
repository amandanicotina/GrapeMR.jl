color_palette(var::AbstractArray) = ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=length(var))))
color_palette(var::Int) = ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=var)))

function plot_config()
    p = plot(
        guidefontsize=12,
        legendfontsize=12,
        tickfontsize=12,
        titlefontsize=15,
        framestyle=:box,
        grid=false,
    )
    return p
end

"""

2D Magnetization Plots

"""
function plot_transverse_magnetization(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    p = plot_config()
    pTrans = plot!(p,
        title="Transverse Magnetization",
        xlabel="Mx",
        ylabel="My")

    max_label_plotted = false
    min_label_plotted = false
    no_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mx = m[2, :]
        My = m[3, :]
        if s.target == "max"
            plot!(pTrans, Mx, My,
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=1.5)
            scatter!(pTrans, [Mx[end]], [My[end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(pTrans, Mx, My,
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=1.5)
            scatter!(pTrans, [Mx[end]], [My[end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(pTrans, Mx, My,
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=1.5)
            scatter!(pTrans, [Mx[end]], [My[end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end
    # xlims!(-1.01, 0.05)
    # ylims!(-0.01, 0.03)

    return pTrans
end


function plot_magnetization_2D(isos::Vector{Isochromat})
    colors = color_palette(10)
    p = plot_config()
    pMag2D = plot!(p,
        title="Magnetization",
        xlabel="Transverse",
        ylabel="Longitudinal")

    max_label_plotted = false
    min_label_plotted = false
    no_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 + m[3, :] .^ 2)
        if s.target == "max"
            plot!(pMag2D, Mxy, m[4, :],
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=1)
            scatter!(pMag2D, [Mxy[end]], [m[4, end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(pMag2D, Mxy, m[4, :],
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=1)
            scatter!(pMag2D, [Mxy[end]], [m[4, end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(pMag2D, Mxy, m[4, :],
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=1)
            scatter!(pMag2D, [Mxy[end]], [m[4, end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end
    xlims!(-0.05, 1.0)
    ylims!(-1.1, 1.1)

    return pMag2D
end


"""

3D Magnetization Plots

"""

function plot_magnetization_3D(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    colors = color_palette(1:4)
    p = plot_config()
    pMag3D = plot!(p,
        title="Magnetization",
        xlabel="Mx",
        ylabel="My",
        zlabel="Mz")
    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        if s.target == "max"
            plot3d!(pMag3D, m[2, :], m[3, :], m[4, :],
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[1], lw=1.5)
            scatter!(pMag3D, [m[2, end]], [m[3, end]], [m[4, end]], label=false, color=colors[1], markersize=4)
            max_label_plotted = true
        elseif s.target == "min"
            plot3d!(pMag3D, m[2, :], m[3, :], m[4, :],
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[3], lw=1.5)
            scatter!(pMag3D, [m[2, end]], [m[3, end]], [m[4, end]], label=false, color=colors[3], markersize=4)
            min_label_plotted = true
        else
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end

    return pMag3D
end



function plot_magnetization_time(iso::Isochromat, t::Float64)
    m = iso.magnetization.dynamics
    s = iso.spin
    time = range(0.0, t, length=length(m[1, :]))

    colors = color_palette(1:10)
    p = plot_config()
    pMagTime = plot!(p,
        title="Magnetization - $(s.label) - $(s.target)",
        xlabel="Transverse",
        ylabel="Longitudinal")

    Mx, My, Mz = m[2, 1:end], m[3, 1:end], m[4, 1:end]
    # Mx
    plot!(pMagTime, time, Mx, label="Mx", color=colors[2], lw=2)
    scatter!(pMagTime, [time[end]], [Mx[end]], label=false, color=colors[2], markersize=5)
    #My
    plot!(pMagTime, time, My, label="My", color=colors[9], lw=2)
    scatter!(pMagTime, [time[end]], [My[end]], label=false, color=colors[9], markersize=5)
    # Mz
    plot!(pMagTime, time, Mz, label="Mz", color=colors[5], lw=2)
    scatter!(pMagTime, [time[end]], [Mz[end]], label=false, color=colors[5], markersize=5)

    return pMagTime
end


function plot_transverse_time(isos::Vector{Isochromat}, t::Float64)
    colors = color_palette(1:10)
    p = plot_config()
    pTransTime = plot!(p,
        title="Transverse Magnetization",
        xlabel="t [sec]",
        ylabel="Mxy")

    max_label_plotted = false
    min_label_plotted = false
    no_label_plotted = false

    mag = isos[1].magnetization.dynamics
    time = range(0.0, t, length=length(mag[1, :]))

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 + m[3, :] .^ 2)
        if s.target == "max"
            plot!(pTransTime, time, Mxy,
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=2)
            scatter!(pTransTime, [time[end]], [Mxy[end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(pTransTime, time, Mxy,
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=2)
            scatter!(pTransTime, [time[end]], [Mxy[end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(pTransTime, time, Mxy,
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=2)
            scatter!(pTransTime, [time[end]], [Mxy[end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end

    return pTransTime
end


function plot_longitudinal_time(isos::Vector{Isochromat}, t::Float64)
    colors = color_palette(1:10)
    p = plot_config()
    pLongTime = plot!(p,
        title="Longitudinal Magnetization",
        xlabel="t [sec]",
        ylabel="Mz")

    max_label_plotted = false
    min_label_plotted = false
    no_label_plotted = false

    mag = isos[1].magnetization.dynamics
    time = range(0.0, t, length=length(mag[1, :]))

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mz = m[4, :]
        if s.target == "max"
            plot!(pLongTime, time, Mz,
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=2)
            scatter!(pLongTime, [time[end]], [Mz[end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(pLongTime, time, Mz,
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=2)
            scatter!(pLongTime, [time[end]], [Mz[end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(pLongTime, time, Mz,
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=2)
            scatter!(pLongTime, [time[end]], [Mz[end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end

    return pLongTime
end





###
# bSSFP stuff
###
function plot_magnetization_target_3D(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    colors = color_palette(1:4)

    # Steady State
    ss = steady_state_matrix(iso)
    Mx_ss, My_ss, Mz_ss = ss.x, ss.y, ss.z

    # Magnetization
    Mx = m[2, :]
    My = m[3, :]
    Mz = m[4, :]

    # Plot
    p = plot3d(title="Magnetization Dynamics - ",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=true)

    plot!(p, Mx, My, Mz, label="$(s.target) - $(s.label)", color=colors[1], lw=2)
    scatter!([Mx[end]], [My[end]], [Mz[end]], label=false, color=colors[1], markersize=5)
    scatter!([Mx_ss], [My_ss], [Mz_ss], label="Steady-state target", color=colors[3], markersize=5)
    xlabel!("Mx")
    ylabel!("My")
    zlabel!("Mz")

    return p
end

function plot_Mtrans_offset_ss(isos::Vector{Isochromat})
    colors = color_palette(1:4)
    label_plotted = false

    p = plot(title="Offset Profile",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 .+ m[3, :] .^ 2)
        b0 = s.B0inho
        # Steady State
        ss = steady_state_matrix(iso)
        Mxy_ss = sqrt((ss.x)^2 + (ss.y)^2)

        scatter!(p, [b0], [Mxy[end]], label=label_plotted ? false : "GrapeMR", color=colors[2], markersize=8)
        scatter!([b0], [Mxy_ss], color=colors[4], label=label_plotted ? false : "Steady State")
        label_plotted = true
    end

    return p
end



# Revisit this target plot functions
function plot_magnetization_target(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    p = plot_config()
    pMagTar = plot!(p,
        title="Magnetization",
        xlabel="Mxy",
        ylabel="Mz")

    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 .+ m[3, :] .^ 2)
        if s.target == "max"
            # Steady State
            ss = steady_state_matrix(iso)
            Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
            plot!(pMagTar, Mxy, m[4, :], label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=2,
                xlims=[-1.0, 1.0],
                ylims=[-1.0, 1.1])
            scatter!(pMagTar, [Mxy[end]], [m[4, end]], label=false, color=colors[2], markersize=4)
            scatter!([Mxy_ss], [Mz_ss], color=colors[9], label=max_label_plotted ? false : "Steady State $(s.target)")
            max_label_plotted = true

        elseif s.target == "min"
            # Steady State
            ss = steady_state_matrix(iso)
            Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
            plot!(pMagTar, Mxy, m[4, :], label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[1], lw=2,
                xlims=[-1.0, 1.0],
                ylims=[-1.0, 1.1])
            scatter!(pMagTar, [Mxy[end]], [m[4, end]], label=false, color=colors[1], markersize=4)
            scatter!([Mxy_ss], [Mz_ss], color=colors[6], label=min_label_plotted ? false : "Steady State $(s.target)")
            min_label_plotted = true
        else
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end

    return pMagTar
end