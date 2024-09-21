"""

2D Magnetization Plots

"""
function plot_transverse_magnetization(isos::Vector{Isochromat})
    pTrans = plot()
    #xlims=[-0.5, 0.5], ylims=[0.5, 1.0])
    for iso ∈ isos
        m = iso.magnetization.dynamics
        Mx = m[2,:]
        My = m[3,:]
        Mt = Mx + im*My
        plot!(pTrans, Mt, color = 2, label = false, lw = 0.5)
        scatter!(pTrans, [Mt[end]], color = 2, label = false)
    end
    return pTrans
end


function plot_magnetization_2D(isos::Vector{Isochromat})
    p = plot(titlefontsize = 12, framestyle=:box)

    max_label_plotted = false
    min_label_plotted = false

    for iso ∈ isos
        m = iso.magnetization.dynamics
        s = iso.spin
        Mxy = sqrt.(m[2,:].^2 + m[3,:].^2)
        if s.target == "max" 
            plot!(p, Mxy, m[4, :], 
            label = max_label_plotted ? false : "$(s.target) - $(s.label)", color = cgrad(:viridis)[128], lw = 1)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = cgrad(:viridis)[128], markersize = 5)
            max_label_plotted = true
            
        elseif s.target == "min" 
            plot!(p, Mxy, m[4, :], 
            label = min_label_plotted ? false : "$(s.target) - $(s.label)", color = :viridis, lw = 1)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color =  :viridis, markersize = 5)
            min_label_plotted = true
        else
            plot!(p, Mxy, m[4, :], 
            label = min_label_plotted ? false : "$(s.target) - $(s.label)", color = "black", lw = 1)
            scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = "black", markersize = 5)
            min_label_plotted = true        
        end
    end
    xlims!(-0.05, 1.0)
    ylims!(-1.1, 1.1)
    xlabel!("Mxy")
    ylabel!("Mz")    
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

function plot_magnetization_targetB0(isos::Vector{Isochromat})
    p = plot()

    for iso ∈ isos
        m   = iso.magnetization.dynamics
        s   = iso.spin
        Mxy = sqrt.(m[2,:].^2 .+ m[3,:].^2)
        # Steady State
        ss = steady_state_matrix(s)
        Mxy_ss, Mz_ss = sqrt((ss.x)^2 + (ss.y)^2), ss.z
        plot!(p, Mxy, m[4, :], label = "$(s.target) - $(s.label)", color = 1, lw = 2,
            xlims = [-1.0, 1.0],
            ylims = [-1.0, 1.1],
            xlabel = "Mxy",
            ylabel = "Mz",
            title = "Magnetization Dynamics",
            titlefontsize = 12)
        scatter!(p, [Mxy[end]], [m[4, end]], label = false, color = 1, markersize = 4)
        scatter!([Mxy_ss], [Mz_ss], color = 3, label = "Steady State $(s.target)")
    end
    return p
end

function plot_Mtrans_offset_ss(isos::Vector{Isochromat})
    label_plotted = false

    p = plot(xlabel = "Offset [Hz]",
        ylabel = "Mxy",
        title = "Offset Profile",
        titlefontsize = 12)

    for iso ∈ isos
        m   = iso.magnetization.dynamics
        s   = iso.spin
        Mxy = sqrt.(m[2,:].^2 .+ m[3,:].^2)
        b0 = s.B0inho
        # Steady State
        ss = steady_state_matrix(s)
        Mxy_ss = sqrt((ss.x)^2 + (ss.y)^2)

        scatter!(p, [b0], [Mxy[end]], label = label_plotted ? false : "GrapeMR", color = 1, markersize = 8)
        scatter!([b0], [Mxy_ss], color = 3, label = label_plotted ? false : "Steady State")
        label_plotted = true
    end

    return p
end

"""

3D Magnetization Plots

"""

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
    #if s.target == "max"
        scatter!([Mx_ss], [My_ss], [Mz_ss], label = "Steady State", color = 3, markersize = 3)
    #end
        xlabel!("Mx")
        ylabel!("My")
        zlabel!("Mz")
        title!("Magnetization Dynamics")

    return p
end
