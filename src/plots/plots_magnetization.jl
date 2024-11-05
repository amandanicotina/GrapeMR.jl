"""
    color_palette(var::Union{AbstractArray, Int})

Generates a color palette using the rainbow color scheme for the specified array of indices or integer count.

# Arguments
- `var::Union{AbstractArray, Int}`: Array or integer specifying the number of colors required.

# Returns
- A `ColorScheme` object with colors matching the length or count of `var`.
"""
function color_palette(var::Union{AbstractArray, Int})
    n = isa(var, Int) ? var : length(var)
    return ColorScheme(get(ColorSchemes.rainbow, range(0.0, 1.0, length=n)))
end

"""
    plot_config()

Sets up a default plot configuration for consistency in appearance across plots.

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
        return (colors[5], "$(s.label)")
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
"""
function plot_transverse_magnetization(isos::Vector{Isochromat})
    colors = color_palette(1:10)
    pTrans = initialize_plot("Transverse Magnetization", "Mx", "My")

    labels_shown = Set()
    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mx, My = m[2, :], m[3, :]
        color, label = get_target_properties(s, colors)
        
        # Plot the trajectory
        plot!(pTrans, Mx, My, label = s.target ∈ labels_shown ? false : label, color = color, lw = 1.5)
        scatter!(pTrans, [Mx[end]], [My[end]], label = false, color = color, markersize = 5)
        push!(labels_shown, s.target)
    end

    return pTrans
end

"""
    plot_magnetization_2D(isos::Vector{Isochromat})

Plots the 2D magnetization, showing transverse (Mxy) vs longitudinal (Mz) magnetization.

# Arguments
- `isos::Vector{Isochromat}`: A vector of `Isochromat` objects containing magnetization and spin data.

# Returns
- `pMag2D`: A plot object displaying the transverse (Mxy) vs longitudinal (Mz) magnetization.
"""
function plot_magnetization_2D(isos::Vector{Isochromat})
    colors = color_palette(10)
    pMag2D = initialize_plot("Magnetization", "Transverse", "Longitudinal")

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

    plot!(pMagTime, time, m[2, :], label = "Mx", color = colors[2], lw = 2)
    scatter!(pMagTime, [time[end]], [m[2, end]], label = false, color = colors[2], markersize = 5)
    plot!(pMagTime, time, m[3, :], label = "My", color = colors[9], lw = 2)
    scatter!(pMagTime, [time[end]], [m[3, end]], label = false, color = colors[9], markersize = 5)
    plot!(pMagTime, time, m[4, :], label = "Mz", color = colors[5], lw = 2)
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

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2, :] .^ 2 + m[3, :] .^ 2)
        color, label = get_target_properties(s, colors)

        plot!(pTransTime, time, Mxy, label = s.target ∈ labels_shown ? false : label, color = color, lw = 2)
        scatter!(pTransTime, [time[end]], [Mxy[end]], label = false, color = color, markersize = 5)
        push!(labels_shown, s.target)
    end

    return pTransTime
end

"""
    plot_longitudinal_time(isos::Vector{Isochromat}, t::Float64)

Plots the longitudinal magnetization (Mz) as a function of time for a set of isochromats.

# Arguments
- `isos::Vector{Isochromat}`: A vector of `Isochromat` objects containing magnetization and spin data.
- `t::Float64`: Total time duration for the plot.

# Returns
- `pLongTime`: A plot object displaying Mz as a function of time for each isochromat.
"""
function plot_longitudinal_time(isos::Vector{Isochromat}, t::Float64)
    colors = color_palette(1:10)
    pLongTime = initialize_plot("Longitudinal Magnetization", "t [sec]", "Mz")

    labels_shown = Set()
    mag = isos[1].magnetization.dynamics
    time = range(0.0, t, length = length(mag[1, :]))

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mz = m[4, :]
        color, label = get_target_properties(s, colors)

        plot!(pLongTime, time, Mz, label = s.target ∈ labels_shown ? false : label, color = color, lw = 2)
        scatter!(pLongTime, [time[end]], [Mz[end]], label = false, color = color, markersize = 5)
        push!(labels_shown, s.target)
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
