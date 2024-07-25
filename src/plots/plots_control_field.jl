"""

Control Field Plots

"""

function plot_control_fields(cf::ControlField)
    time = range(0.0, cf.t_control, length = length(cf.B1x))
    Bx = cf.B1x
    By = cf.B1y

    B1 = vec(Bx + im * By)

    p_Bx = plot(time, abs.(B1), linewidth=2, label=false, ylabel="|B1| [Hz]", title="Control Fields", titlefontsize=15)
    p_By = plot(time, angle.(B1), linewidth=2, label=false, ylabel="ϕ [rad]", xlabel="t [s]")
    xticks_values = [-π, -π/2, 0, π/2, π]
    xticks_labels = ["-π", "-π/2", "0", "π/2", "π"]
    p_By = plot!(p_By, xticks=(xticks_values, xticks_labels))

    p = plot(p_Bx, p_By, layout = (2,1))

    return p
end


function plot_control_fields_tesla(cf::ControlField)
    time = range(0.0, cf.t_control, length = length(cf.B1x))
    Bx = cf.B1x./γ_¹H
    By = cf.B1y./γ_¹H

    B1 = vec(1e6*(Bx + im*By))
    
    p_Bx = plot(time, abs.(B1), linewidth=2, label=false, ylabel="|B1| [μT]", title="Control Fields", titlefontsize=15)
    p_By = plot(time, angle.(B1), linewidth=2, label=false, ylabel="ϕ [rad]", xlabel="t [s]")
    xticks_values = [-π, -π/2, 0, π/2, π]
    xticks_labels = ["-π", "-π/2", "0", "π/2", "π"]
    p_By = plot!(p_By, xticks=(xticks_values, xticks_labels))

    p = plot(p_Bx, p_By, layout=(2,  1))

    return p 
end

function plot_control_fields_phase(cf::ControlField; ψ::Float64 = π)
    ux = cf.B1x
    uy = cf.B1y
    Ux = vec(ux*cos(ψ) .+ uy*sin(ψ)) 
    Uy = vec(uy*cos(ψ) .- ux*sin(ψ))

    tc = cf.t_control
    t = range(0.0, tc, length = length(ux))
    p_Bx = plot(t, Ux, linewidth=2, label=false, ylabel="B1x", title="Control Fields", titlefontsize=15)
    p_By = plot(t, Uy, linewidth=2, label=false, ylabel="B1y", xlabel="t [s]")
    # xticks_values = [-π, -π/2, 0, π/2, π]
    # xticks_labels = ["-π", "-π/2", "0", "π/2", "π"]
    # p_By = plot!(p_By, xticks=(xticks_values, xticks_labels))
    p = plot(p_Bx, p_By, layout = (2,1))
    return p
end


function plot_magnetization_control_field(cf::ControlField, isos::Vector{Isochromat})
    max_label_plotted = false
    min_label_plotted = false
    pTrans = plot(title = "Transverse Magnetization",
            titlefontsize = 12,
            xlabel="Mx", ylabel="My",
            framestyle=:box)
    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        if s.target == "max" 
            Mt =  m[2,:] + im*m[3,:]
            Mx = real.(Mt)
            My = imag.(Mt)
            plot!(pTrans, Mx, My, color = cgrad(:viridis)[128], lw = 1.5,
            label = max_label_plotted ? false : "$(s.target) - $(s.label)")
            scatter!(pTrans, [Mx[end]], [My[end]], label = false, color = cgrad(:viridis)[128], markersize=4)
            max_label_plotted = true

        elseif s.target == "min" 
            Mt =  m[2,:] + im*m[3,:]
            Mx = real.(Mt)
            My = imag.(Mt)
            plot!(pTrans, Mx, My, color = "black", lw = 1.5,
            label = min_label_plotted ? false : "$(s.target) - $(s.label)") 
            scatter!(pTrans, [Mx[end]], [My[end]], label = false, color = "black", markersize=4)
            min_label_plotted = true
        end
    end

    time = range(0.0, cf.t_control, length = length(cf.B1x))*1e3
    Bx = cf.B1x
    By = cf.B1y
    B1 = vec(Bx + im * By)

    p_Bx = plot(time, abs.(B1), linewidth=2, grid = false, label=false, color=:viridis, framestyle=:box, ylabel="|B1| [Hz]", title="Control Fields", titlefontsize=12)
    p_By = plot(time, angle.(B1), linewidth=2, grid = false, label=false, color=:viridis, framestyle=:box, ylabel="ϕ [rad]", xlabel="t [ms]")

    l = @layout [[a{0.7h}; b{0.3h}] c]
    pMagCf = plot(p_Bx, p_By, pTrans, layout = l)
    display(pMagCf)

end



