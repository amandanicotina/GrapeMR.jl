"""
    plot_control_fields(cf::ControlField; unit::String = "Hz")

Plot the control fields of an RF pulse in different units.

# Arguments
- `cf::ControlField`: The control field object containing the RF waveform data (`B1x`, `B1y` components).
- `unit::String`: The unit for the plot. Supported units are:
    - `"Hz"`: Plots `B1x` and `B1y` in Hertz.
    - `"rad/s"`: Converts the control fields to rads/sec and plots the amplitude and phase.
    - `"Tesla"`: Converts the control fields to Tesla using the gyromagnetic ratio for ¹H and plots in microtesla (µT).

# Returns
- A plot object displaying the control fields in the specified units with amplitude and phase information.
"""
function plot_control_fields(cf::ControlField; unit::String = "Hz")
    colors = color_palette(10)
    time = range(0.0, cf.t_control, length = length(cf.B1x))
    Bx, By = vec(cf.B1x), vec(cf.B1y)
    ylabel_B1, ylabel_phase = "", ""

    # Adjust fields based on unit
    if unit == "Hz"
        ylabel_B1, ylabel_phase = "Bx [Hz]", "By [Hz]"
    elseif unit == "rad/s"
        Bx, By = 2π .* Bx, 2π .* By  # Convert Hz to rad/s
        B1 = vec(Bx + im * By)
        Bx, By = abs.(B1), angle.(B1)
        ylabel_B1, ylabel_phase = "|B1| [rad/s]", "ϕ [rad]"
    elseif unit == "Tesla"
        γ_¹H = 42.577e6  # Gyromagnetic ratio in Hz/T for ¹H
        Bx, By = Bx ./ γ_¹H, By ./ γ_¹H
        B1 = vec(Bx + im * By)
        Bx, By = abs.(B1) * 1e6, angle.(B1)
        ylabel_B1, ylabel_phase = "|B1| [μT]", "ϕ [rad]"
    else
        error("Invalid unit. Supported units are 'Hz', 'rad/s', and 'Tesla'.")
    end

    # Use initialize_plot to set up plots for each component
    p_Bx = initialize_plot("Control Fields", "t [s]", ylabel_B1)
    plot!(p_Bx, time, Bx, color = colors[3], lw = 2, label = false)

    p_By = initialize_plot("", "t [s]", ylabel_phase)
    plot!(p_By, time, By, color = colors[3], lw = 2, label = false)

    # Combine plots into a single layout
    p = plot(p_Bx, p_By, layout = (2, 1))
    return p
end


"""
    plot_control_fields_phase_shift(cf::ControlField; ψ::Float64 = π)

Plot the control fields of an RF pulse after applying a phase shift.

# Arguments
- `cf::ControlField`: The control field object containing the RF waveform data (`B1x`, `B1y` components).
- `ψ::Float64`: The phase shift to apply to the control fields (default is π radians).

# Returns
- A plot object displaying the control fields with the specified phase shift.
"""
function plot_control_fields_phase_shift(cf::ControlField; ψ::Real = π)
    colors = color_palette(10)
    time = range(0.0, cf.t_control, length = length(cf.B1x))

    # Apply phase shift to control fields
    Ux = vec(cf.B1x * cos(ψ) .+ cf.B1y * sin(ψ))
    Uy = vec(cf.B1y * cos(ψ) .- cf.B1x * sin(ψ))

    # Initialize plots using helper function
    p_Bx = initialize_plot("Control Fields with Phase Shift", "t [s]", "B1x")
    plot!(p_Bx, time, Ux, color = colors[3], lw = 2, label = false)

    p_By = initialize_plot("", "t [s]", "B1y")
    plot!(p_By, time, Uy, color = colors[3], lw = 2, label = false)

    # Combine plots into a single layout
    p = plot(p_Bx, p_By, layout = (2, 1))
    return p
end


"""
    plot_magnetization_control_field(cf::ControlField, isos::Vector{Isochromat})

Plot the magnetization trajectory for a set of isochromats and the corresponding control field.

# Arguments
- `cf::ControlField`: The control field object containing the RF waveform data (`B1x`, `B1y` components).
- `isos::Vector{Isochromat}`: A vector of isochromat objects representing different spin systems.

# Returns
- A plot object displaying:
    1. The magnetization trajectory in the transverse plane.
    2. The control field amplitude and phase over time.
"""
function plot_magnetization_control_field(cf::ControlField, isos::Vector{<:Isochromat})
    colors = color_palette(10)
    labels_shown = Set()

    # Initialize magnetization trajectory plot
    pMag = initialize_plot("Magnetization Trajectory", "Transverse", "Longitudinal")
    xlims!(-0.02, 1.01)
    ylims!(-1.01, 1.01)

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy, Mz = sqrt.(m[2, :] .^ 2 + m[3, :] .^ 2), m[4, :]
        color, label, α = get_target_properties(s, colors)

        # Plot trajectory and final point
        plot!(pMag, Mxy, Mz, label = s.target ∈ labels_shown ? false : label, color = color, lw = 4, alpha = α)
        scatter!(pMag, [Mxy[end]], [Mz[end]], label = false, color = color, markersize = 8)
        push!(labels_shown, s.target)
    end

    # Prepare time vector and control field components
    time = range(0.0, cf.t_control, length = length(cf.B1x)) * 1e3  # Convert to milliseconds

    # Initialize control field plots using helper function
    p_Bx = initialize_plot("Control Fields", "t [ms]", "Bx [Hz]")
    plot!(p_Bx, time, vec(cf.B1x), label = false, color = colors[5], lw = 4)

    p_By = initialize_plot("", "t [ms]", "By [Hz]")
    plot!(p_By, time, vec(cf.B1y), label = false, color = colors[5], lw = 4)

    # Arrange plots into the final layout
    l = @layout [[a{0.5h}; b{0.5h}] c]
    pMagCf = plot(p_Bx, p_By, pMag, layout = l, size = (800, 600))

    return pMagCf
end


"""
    plot_cost_values(cost::Vector{Float64}, gp::GrapeParams)

Plot the cost function convergence over iterations during the GRAPE optimization process.

# Arguments
- `cost::Vector{Float64}`: A vector containing the cost values at each iteration.
- `gp::GrapeParams`: The parameters of the GRAPE optimization, including the cost function name.

# Returns
- A plot object displaying the convergence of the cost function over the iterations.
"""
function plot_cost_values(cost::Vector{Float64}, gp::GrapeParams)
    colors = color_palette(10)
    cmin = round(cost[end], digits = 3)
    iter = range(0, stop = length(cost), length = length(cost))
    
    # Initialize the plot using the helper function
    p = initialize_plot("Cost Function Convergence", "Iterations", "Cost Value")
    
    # Plot the cost values over iterations
    plot!(p, iter, cost, label = " " * string(gp.cost_function), lw = 3, color = colors[2])
    
    # Highlight the final cost value with a scatter point
    scatter!(p, [iter[end]], [cost[end]], color = colors[2], markersize = 7, label = "Final Cost = $cmin")

    return p
end
