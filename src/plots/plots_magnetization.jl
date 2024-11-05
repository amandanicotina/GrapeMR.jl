color_palette(var::AbstractArray) = ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=length(var))))
color_palette(var::Int) = ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=var)))

"""

2D Magnetization Plots

<<<<<<< Updated upstream
=======
# Returns
- `p`: A plot object with font sizes, frame style, and grid style preconfigured.
"""
function plot_config()
    plot(
        guidefontsize = 15,
        legendfontsize = 15,
        tickfontsize = 10,
        titlefontsize = 17,
        framestyle = :box,
        grid = false,
    )
end

"""
    initialize_plot(title, xlabel, ylabel; zlabel=nothing)

Creates a plot object with common settings and title and axis labels.

# Arguments
- `title::String`: Title of the plot.
- `xlabel::String`: Label for the x-axis.
- `ylabel::String`: Label for the y-axis.
- `zlabel::Union{String, Nothing}`: Label for the z-axis (optional, only for 3D plots).

# Returns
- `p`: A plot object with preconfigured labels and settings.
"""
function initialize_plot(title::String, xlabel::String, ylabel::String; zlabel::Union{String, Nothing} = nothing)
    p = plot_config()
    plot!(p, title = title, xlabel = xlabel, ylabel = ylabel)
    if !isnothing(zlabel)
        zlabel!(zlabel)
    end
    return p
end

"""
    get_target_properties(s, colors)

Determines the color and label for a given spin target ("max", "min", or default).

# Arguments
- `s::Spin`: A spin object with target and label properties.
- `colors::Vector`: A vector of colors from the color palette.

# Returns
- A tuple `(color, label)` with the color and label for the given spin target.
"""
function get_target_properties(s, colors)
    if s.target == "max"
        return (colors[2], "$(s.target) - $(s.label)")
    elseif s.target == "min"
        return (colors[9], "$(s.target) - $(s.label)")
    else
        return (colors[1], "$(s.label)")
    end
end

# 2D and 3D Magnetization Plotting Functions

"""
    plot_transverse_magnetization(isos::Vector{Isochromat})

Plots the transverse magnetization (Mx and My) for a set of isochromats, with colors indicating different spin targets.

# Arguments
- `isos::Vector{Isochromat}`: A vector of `Isochromat` objects containing magnetization and spin information.

# Returns
- `pTrans`: A plot object displaying the transverse magnetization (Mx vs My) for each isochromat.
>>>>>>> Stashed changes
"""
function plot_transverse_magnetization(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    pTrans = plot(title="Transverse Magnetization",
        xlabel="Mx",
        ylabel="My",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

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
    p = plot(title="Magnetization",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

<<<<<<< Updated upstream
    max_label_plotted = false
    min_label_plotted = false
    no_label_plotted = false
=======
    labels_shown = Set()
    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy, Mz = sqrt.(m[2, :] .^ 2 + m[3, :] .^ 2), m[4, :]
        color, label = get_target_properties(s, colors)
        
        plot!(pMag2D, Mxy, Mz, label = s.target ∈ labels_shown ? false : label, color = color, lw = 1)
        scatter!(pMag2D, [Mxy[end]], [Mz[end]], label = false, color = color, markersize = 5)
        push!(labels_shown, s.target)
    end

    xlims!(pMag2D, -0.05, 1.0)
    ylims!(pMag2D, -1.1, 1.1)

    return pMag2D
end

"""
    plot_magnetization_3D(isos::Vector{Isochromat})

Plots the 3D magnetization (Mx, My, Mz) for a set of isochromats.

# Arguments
- `isos::Vector{Isochromat}`: A vector of `Isochromat` objects containing magnetization and spin data.

# Returns
- `pMag3D`: A plot object displaying the 3D magnetization components (Mx, My, Mz) for each isochromat.
"""
function plot_magnetization_3D(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    pMag3D = initialize_plot("Magnetization", "Mx", "My", zlabel = "Mz")

    labels_shown = Set()
    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        color, label = get_target_properties(s, colors)
        
        plot3d!(pMag3D, m[2, :], m[3, :], m[4, :], label = s.target ∈ labels_shown ? false : label, color = color, lw = 1.5)
        scatter!(pMag3D, [m[2, end]], [m[3, end]], [m[4, end]], label = false, color = color, markersize = 4)
        push!(labels_shown, s.target)
    end

    return pMag3D
end

"""
    plot_magnetization_time(iso::Isochromat, t::Float64)

Plots the time evolution of magnetization components (Mx, My, Mz) for a single isochromat.

# Arguments
- `iso::Isochromat`: An `Isochromat` object containing magnetization dynamics.
- `t::Float64`: Total time duration of the plot.

# Returns
- `pMagTime`: A plot object displaying Mx, My, and Mz as functions of time.
"""
function plot_magnetization_time(iso::Isochromat, t::Float64)
    colors = color_palette(1:10)
    pMagTime = initialize_plot("Magnetization - $(iso.spin.label) - $(iso.spin.target)", "Time [s]", "Magnetization")

    m = iso.magnetization.dynamics
    time = range(0.0, t, length = length(m[1, :]))

    plot!(pMagTime, time, m[2, :], label = "Mx", color = colors[1], lw = 2)
    scatter!(pMagTime, [time[end]], [m[2, end]], label = false, color = colors[2], markersize = 5)
    plot!(pMagTime, time, m[3, :], label = "My", color = colors[5], lw = 2)
    scatter!(pMagTime, [time[end]], [m[3, end]], label = false, color = colors[9], markersize = 5)
    plot!(pMagTime, time, m[4, :], label = "Mz", color = colors[9], lw = 2)
    scatter!(pMagTime, [time[end]], [m[4, end]], label = false, color = colors[5], markersize = 5)

    return pMagTime
end

"""
    plot_transverse_time(isos::Vector{Isochromat}, t::Float64)

Plots the transverse magnetization (Mxy) as a function of time for a set of isochromats.

# Arguments
- `isos::Vector{Isochromat}`: A vector of `Isochromat` objects containing magnetization and spin data.
- `t::Float64`: Total time duration for the plot.

# Returns
- `pTransTime`: A plot object displaying Mxy as a function of time for each isochromat.
"""
function plot_transverse_time(isos::Vector{Isochromat}, t::Float64)
    colors = color_palette(1:10)
    pTransTime = initialize_plot("Transverse Magnetization", "t [sec]", "Mxy")

    labels_shown = Set()
    mag = isos[1].magnetization.dynamics
    time = range(0.0, t, length = length(mag[1, :]))
>>>>>>> Stashed changes

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 + m[3, :] .^ 2)
        if s.target == "max"
            plot!(p, Mxy, m[4, :],
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=1)
            scatter!(p, [Mxy[end]], [m[4, end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(p, Mxy, m[4, :],
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=1)
            scatter!(p, [Mxy[end]], [m[4, end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(p, Mxy, m[4, :],
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=1)
            scatter!(p, [Mxy[end]], [m[4, end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end
    xlims!(-0.05, 1.0)
    ylims!(-1.1, 1.1)
    xlabel!("Transverse")
    ylabel!("Longitudinal")

    return p
end

# Revisit this target plot functions
function plot_magnetization_target(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    p = plot(title="Magnetization Dynamics",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

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
            plot!(p, Mxy, m[4, :], label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=2,
                xlims=[-1.0, 1.0],
                ylims=[-1.0, 1.1],
                xlabel="Mxy",
                ylabel="Mz")
            scatter!(p, [Mxy[end]], [m[4, end]], label=false, color=colors[2], markersize=4)
            scatter!([Mxy_ss], [Mz_ss], color=colors[9], label=max_label_plotted ? false : "Steady State $(s.target)")
            max_label_plotted = true

        elseif s.target == "min"
            # Steady State
            ss = steady_state_matrix(iso)
            Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
            plot!(p, Mxy, m[4, :], label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[1], lw=2,
                xlims=[-1.0, 1.0],
                ylims=[-1.0, 1.1],
                xlabel="Mxy",
                ylabel="Mz")
            scatter!(p, [Mxy[end]], [m[4, end]], label=false, color=colors[1], markersize=4)
            scatter!([Mxy_ss], [Mz_ss], color=colors[6], label=min_label_plotted ? false : "Steady State $(s.target)")
            min_label_plotted = true
        else
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end

    return p
end

function plot_magnetization_targetB0(isos::Vector{Isochromat})
    colors = color_palette(Int(2 * length(isos)))
    p = plot(title="Magnetization Dynamics - ",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

    for (i, iso) ∈ enumerate(isos)
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 .+ m[3, :] .^ 2)
        # Steady State
        ss = steady_state_matrix(iso)
        Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
        plot!(p, Mxy, m[4, :], label="ΔB0 = $(s.B0inho)Hz", color=colors[i], lw=2,
            xlims=[-1.0, 1.0],
            ylims=[-1.0, 1.1],
            xlabel="Mxy",
            ylabel="Mz",
            title="Magnetization Dynamics",
            titlefontsize=12)
        scatter!(p, [Mxy[end]], [m[4, end]], label=false, color=colors[i], markersize=4)
        scatter!([Mxy_ss], [Mz_ss], color=colors[2*i], label="Steady-state target")
    end
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

"""

3D Magnetization Plots

"""

function plot_magnetization_3D(isos::Vector{Isochromat})
    p = plot3d(title="Magnetization Dynamics - ",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)
    colors = color_palette(1:4)
    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        if s.target == "max"
            plot3d!(p, m[2, :], m[3, :], m[4, :],
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[1], lw=1.5)
            scatter!(p, [m[2, end]], [m[3, end]], [m[4, end]], label=false, color=colors[1], markersize=4)
            max_label_plotted = true
        elseif s.target == "min"
            plot3d!(p, m[2, :], m[3, :], m[4, :],
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[3], lw=1.5)
            scatter!(p, [m[2, end]], [m[3, end]], [m[4, end]], label=false, color=colors[3], markersize=4)
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


function plot_magnetization_time(iso::Isochromat, t::Float64)
    mag = iso.magnetization.dynamics
    spin = iso.spin
    time = range(0.0, t, length=length(mag[1, :]))

    # Create a plot of the magnetization
    p = plot(time, mag[2:end, :]', label=["Mx" "My" "Mz"], lw=2,
        xlabel="t [sec]",
        ylabel="Magnitude",
        title="Magnetization Dynamics - Sample = $(spin.label), Target = $(spin.target)",
        titlefontsize=12,
    )
    return p
end


function plot_transverse_time(isos::Vector{Isochromat}, t::Float64)
    colors = color_palette(1:10)
    pTrans = plot(title="Transverse Magnetization",
        xlabel="t [sec]",
        ylabel="Mxy",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

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
            plot!(pTrans, time, Mxy,
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=2)
            scatter!(pTrans, [time[end]], [Mxy[end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(pTrans, time, Mxy,
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=2)
            scatter!(pTrans, [time[end]], [Mxy[end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(pTrans, time, Mxy,
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=2)
            scatter!(pTrans, [time[end]], [Mxy[end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end

    return pTrans
end


function plot_longitudinal_time(isos::Vector{Isochromat}, t::Float64)
    colors = color_palette(1:10)
    pLong = plot(title="Longitudinal Magnetization",
        xlabel="t [sec]",
        ylabel="Mz",
        guidefontsize=12,
        legendfontsize=10,
        tickfontsize=10,
        titlefontsize=12,
        framestyle=:box,
        grid=false)

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
            plot!(pLong, time, Mz,
                label=max_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[2], lw=2)
            scatter!(pLong, [time[end]], [Mz[end]], label=false, color=colors[2], markersize=5)
            max_label_plotted = true

        elseif s.target == "min"
            plot!(pLong, time, Mz,
                label=min_label_plotted ? false : "$(s.target) - $(s.label)", color=colors[9], lw=2)
            scatter!(pLong, [time[end]], [Mz[end]], label=false, color=colors[9], markersize=5)
            min_label_plotted = true
        else
            plot!(pLong, time, Mz,
                label=min_label_plotted ? false : "$(s.label)", color=colors[5], lw=2)
            scatter!(pLong, [time[end]], [Mz[end]], label=false, color=colors[5], markersize=5)
            no_label_plotted = true
        end
    end

    return pLong
end
