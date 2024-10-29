"""
    plot_control_fields(cf::ControlField; unit::String = "Hz")

Plot the control fields of an RF pulse in different units.

    ### Parameters:
    - `cf::ControlField`: The control field object containing the RF waveform data (`B1x`, `B1y` components).
    - `unit::String`: The unit for the plot. Supported units are:
        - `"Hz"`: Plots `B1x` and `B1y` in Hertz.
        - `"rad/s"`: Converts the control fields to rads/sec and plots the amplitude and phase.
        - `"Tesla"`: Converts the control fields to Tesla using the gyromagnetic ratio for ¹H and plots in microtesla (µT).

    ### Returns:
    - A plot object displaying the control fields in the specified units with amplitude and phase information.
"""
function plot_control_fields(cf::ControlField; unit::String = "Hz")
    colors = color_palette(10)
    time = range(0.0, cf.t_control, length = length(cf.B1x))
    Bx, By = cf.B1x, cf.B1y
    ylabel_B1, ylabel_phase = "", ""

    if unit == "Hz"
        Bx, By = vec(Bx), vec(By)
        ylabel_B1 = "Bx [Hz]"
        ylabel_phase = "By [Hz]"
    elseif unit == "rad/s"
        Bx, By = 2π .* Bx, 2π .* By # Convert Hz to rad/s
        B1 = vec(Bx + im * By)
        Bx, By = abs.(B1), angle.(B1)
        ylabel_B1 = "|B1| [rad/s]"
        ylabel_phase = "ϕ [rad]"
    elseif unit == "Tesla"
        γ_¹H = 42.577e6  # Gyromagnetic ratio in Hz/T for ¹H
        Bx, By = Bx ./ γ_¹H, By ./ γ_¹H
        B1 = vec(Bx + im * By)
        Bx, By = abs.(B1)*1e6, angle.(B1)
        ylabel_B1 = "|B1| [μT]"
        ylabel_phase = "ϕ [rad]"
    else
        error("Invalid unit. Supported units are 'Hz', 'rad/s', and 'Tesla'.")
    end


    p_Bx = plot(time, Bx, color = colors[3], lw=2, label=false, ylabel=ylabel_B1, title="Control Fields", titlefontsize=15)
    p_By = plot(time, By, color = colors[3], lw=2, label=false, ylabel=ylabel_phase, xlabel="t [s]")


    p = plot(p_Bx, p_By, layout=(2, 1),             
            guidefontsize = 12, 
            legendfontsize = 10,
            tickfontsize = 10,
            titlefontsize = 14,
            framestyle = :box, 
            grid = false)

    return p
end


"""
    plot_control_fields_phase_shift(cf::ControlField; ψ::Float64 = π)

Plot the control fields of an RF pulse after applying a phase shift.

    ### Parameters:
    - `cf::ControlField`: The control field object containing the RF waveform data (`B1x`, `B1y` components).
    - `ψ::Float64`: The phase shift to apply to the control fields (default is π radians).

    ### Returns:
    - A plot object displaying the control fields with the specified phase shift.
"""

function plot_control_fields_phase_shift(cf::ControlField; ψ::Float64 = π)
    colors = color_palette(10)
    ux = cf.B1x
    uy = cf.B1y
    Ux = vec(ux*cos(ψ) .+ uy*sin(ψ)) 
    Uy = vec(uy*cos(ψ) .- ux*sin(ψ))

    tc = cf.t_control
    t = range(0.0, tc, length = length(ux))
    p_Bx = plot(t, Ux, color = colors[3], lw=2, label=false, ylabel="B1x", title="Control Fields", titlefontsize=15)
    p_By = plot(t, Uy, color = colors[3], lw=2, label=false, ylabel="B1y", xlabel="t [s]")
    # xticks_values = [-π, -π/2, 0, π/2, π]
    # xticks_labels = ["-π", "-π/2", "0", "π/2", "π"]
    # p_By = plot!(p_By, xticks=(xticks_values, xticks_labels))
    
    p = plot(p_Bx, p_By, layout=(2, 1),             
            guidefontsize = 12, 
            legendfontsize = 10,
            tickfontsize = 10,
            titlefontsize = 14,
            framestyle = :box, 
            grid = false)

    return p
end



"""
    plot_magnetization_control_field(cf::ControlField, isos::Vector{Isochromat})

Plot the magnetization trajectory for a set of isochromats and the corresponding control field.

    ### Parameters:
    - `cf::ControlField`: The control field object containing the RF waveform data (`B1x`, `B1y` components).
    - `isos::Vector{Isochromat}`: A vector of isochromat objects representing different spin systems.

    ### Returns:
    - A plot object displaying:
        1. The magnetization trajectory in the transverse plane.
        2. The control field amplitude and phase over time.
"""
function plot_magnetization_control_field(cf::ControlField, isos::Vector{Isochromat})
    colors = color_palette(10)
    max_label_plotted = false
    min_label_plotted = false
    pMag = plot(title = "Magnetization Trajectory",
            titlefontsize = 12,
            xlims = (-0.05, 1.0),
            ylims = (-1.1, 1.1),
            xlabel = "Transverse",
            ylabel = "Longitudinal",   
            framestyle=:box)

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2,:].^2 + m[3,:].^2)
        if s.target == "max" 
            plot!(pMag, Mxy, m[4, :], 
            label = max_label_plotted ? false : "$(s.target) - $(s.label)", color = colors[2], lw = 2)
            scatter!(pMag, [Mxy[end]], [m[4, end]], label = false, color = colors[2], markersize = 4)
            max_label_plotted = true
            
        elseif s.target == "min" 
            plot!(pMag, Mxy, m[4, :], 
            label = min_label_plotted ? false : "$(s.target) - $(s.label)", color = colors[9], lw = 2)
            scatter!(pMag, [Mxy[end]], [m[4, end]], label = false, color = colors[9], markersize = 4)
            min_label_plotted = true
        else
            plot!(pMag, Mxy, m[4, :], 
            label = min_label_plotted ? false : "$(s.label)", color = colors[2], lw = 2)
            scatter!(pMag, [Mxy[end]], [m[4, end]], label = false, color = colors[2], markersize = 5)
            min_label_plotted = true        
        end
    end

    time = range(0.0, cf.t_control, length = length(cf.B1x))*1e3
    Bx = vec(cf.B1x)
    By = vec(cf.B1y)
    B1 = vec(Bx + im * By)

    p_Bx = plot(time, Bx, color = colors[5], linewidth=2, grid = false, label=false, framestyle=:box, ylabel="Bx [Hz]", title="Control Fields", titlefontsize=12)
    p_By = plot(time, By, color = colors[5], linewidth=2, grid = false, label=false, framestyle=:box, ylabel="By [Hz]", xlabel="t [ms]") #ϕ [rad]

    l = @layout [[a{0.5h}; b{0.5h}] c]
    
    pMagCf = plot(p_Bx, p_By, pMag, layout = l,
                guidefontsize = 12, 
                legendfontsize = 9,
                tickfontsize = 9,
                titlefontsize = 12,
                framestyle = :box, 
                grid = false)

    display(pMagCf)
end



"""
    plot_cost_values(cost::Vector{Float64}, gp::GrapeParams)

Plot the cost function convergence over iterations during the GRAPE optimization process.

    ### Parameters:
    - `cost::Vector{Float64}`: A vector containing the cost values at each iteration.
    - `gp::GrapeParams`: The parameters of the GRAPE optimization, including the cost function name.

    ### Returns:
    - A plot object displaying the convergence of the cost function over the iterations.
"""
function plot_cost_values(cost::Vector{Float64}, gp::GrapeParams)
    colors = color_palette(10)
    cmin = round(cost[end], digits = 3)
    iter = range(0.0, stop = length(cost), length = length(cost))
    p = plot(iter, cost, label = string(gp.cost_function), lw = 2, color = colors[2],
           # ylims = (0.0, 1.0),
            xlabel = "Iterations",
            ylabel = "Cost Value",
            title  = "Cost Function Convergence",
            guidefontsize = 12, 
            legendfontsize = 10,
            tickfontsize = 10,
            titlefontsize = 12,
            framestyle = :box, 
            grid = false)
        scatter!([iter[end]], [cost[end]], color = colors[2],  markersize = 5, label = "Final Cost = $cmin")
        

    return p
end

