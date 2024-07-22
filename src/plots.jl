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
 
function plot_transverse_magnetization(isos::Vector{Isochromat})
    pTrans = plot()
    #xlims=[-0.5, 0.5], ylims=[0.5, 1.0])
    for iso ∈ isos
        m = iso.magnetization.dynamics
        Mx = m[2,:]
        My = m[3,:]
        Mt = Mx + im*My
        plot!(pTrans, Mt, color = 2, label = false)
        scatter!(pTrans, [Mt[end]], color = 2, label = false)
    end
    return pTrans
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
            error(" $(s.target) is not a matching target. Valid targets are max or min")
        end
    end
    xlims!(-0.05, 1.0)
    ylims!(-1.1, 1.1)
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
    #if s.target == "max"
        scatter!([Mx_ss], [My_ss], [Mz_ss], label = "Steady State", color = 3, markersize = 3)
    #end
        xlabel!("Mx")
        ylabel!("My")
        zlabel!("Mz")
        title!("Magnetization Dynamics")

    return p
end

"""

Cost Function plots

"""

function plot_cost_values(cost::Vector{Float64}, gp::GrapeParams)
    p = plot(cost, label = string(gp.cost_function), lw = 2,
    xlabel = "Iterations",
    ylabel = "Cost Value",
    title  = "Cost Function Convergence",
    titlefontsize = 12,
    )

    return p
end


function plot_cost_offset(isos::Vector{Isochromat}, cost::Symbol)
    # Cost function values for all isochromats
    c = cost_function.(isos, cost)

    # Offset frequencies
    ν_ini = isos[1].spin.B0inho
    ν_end = isos[end].spin.B0inho
    ν_len = Int(ceil(length(isos)/2))
    ν = range(ν_ini, stop=ν_end, length=ν_len)

    p = plot(xlabel = "Offset [Hz]",
        ylabel = "Cost Value",
        title  = "Cost Function Offset profile",
        titlefontsize = 12)
        plot!(p, ν, c[1:ν_len], label =  "min", lw = 2)
        plot!(p, ν, c[ν_len+1:end], label = "max", lw = 2)

    return p
end

function plot_cost_offset(spins::Vector{GrapeMR.Spin}, cost::Symbol)
    # Calculate dynamics for new offset range with optimized field
    isos = Vector{Isochromat}()
    for spin ∈ spins
        mag = forward_propagation(grape_output.control_field, spin)
        dyn = GrapeMR.Magnetization(mag)
        iso = Isochromat(dyn, spin)
        push!(isos, iso)
    end
    # Cost function values for all isochromats
    c = cost_function.(isos, cost)

    # Offset frequencies
    ν_ini = isos[1].spin.B0inho
    ν_end = isos[end].spin.B0inho
    ν_len = ceil(length(isos))
    ν = range(ν_ini, stop=ν_end, length=ν_len)

    p = plot(xlabel = "Offset [Hz]",
        ylabel = "Cost Value",
        title  = "Cost Function Offset profile",
        titlefontsize = 12)
        plot!(p, ν, c, label =  "min", lw = 2)
        # plot!(p, ν, c[ν_len+1:end], label = "max", lw = 2)

    return p
end

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


function bohb_params(bohb)
    cost  = []
    t_c   = []
    start = []
    deg   = []
    max   = []
    for i ∈ eachindex(bohb.results)
        push!(cost, bohb.results[i][1])
        push!(t_c, bohb.history[i][1])
        push!(start, bohb.history[i][2])
        push!(deg, bohb.history[i][3])
        push!(max, bohb.history[i][4])
    end
    function get_eps(start, deg, max_iter)
        return start / (1 - (max_iter/2 - 1) / max_iter)^deg
    end
    eps = [get_eps(start, deg, it) for (start, deg, it) in zip(start, deg, max)]
    t_min = bohb.minimizer[1]
    st_min = bohb.minimizer[2]
    d_min = bohb.minimizer[3]
    m_min = bohb.minimizer[4]
    c_min = round(bohb.minimum[1], digits=4)
    order = collect(1:length(t_c))
    p_cost =  scatter(t_c, log.(max), zcolor=cost, 
         markerstrokecolor = :auto, label = false,
         xlabel = "Control time [s]", ylabel = "Iter", colorbar_title="Cost Value",
         title = "Optimization",
         color = :viridis)
         scatter!([t_min], [log.(m_min)], label = "Minimum = $c_min", 
         marker = :star5, markersize = 8, color = :red),
    p_order =  scatter(cost, log.(max), zcolor=order, 
         markerstrokecolor = :auto, label = false, 
         xlabel = "Cost Function", ylabel = "Resources - log scale", colorbar_title="Iterations",
         title = "Hyperparameter Tuning - Random Sampler", 
         color = :viridis),
         scatter!([c_min], [log(m_min)], label = "Minimum = $c_min", 
         marker = :star5, markersize = 8, color = :red)

    return p_cost, p_order
end

# p1, p2 = bohb_params(rand_hopt)

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


function plot_evaluations(bohb)
    plots = []
    for (i, i_param) in enumerate(bohb.params)
        for (j, j_param) in enumerate(bohb.params)
            if i == j
                hist = histogram([val[i] for val in bohb.history], bins=10, xlabel=i_param, ylabel="Number of samples", legend=false)
                push!(plots, hist)
            elseif j > i
                # we don't care about these
                continue
            else
                # plot x,y with colour equivalent to cost value
                xs = [val[j] for val in bohb.history]
                ys = [val[i] for val in bohb.history]
                xy = scatter(
                    xs, ys, zcolor=bohb.results, xlabel=j_param, ylabel=i_param, legend=false,
                    color = :viridis
                )
                scatter!([bohb.minimizer[j]], [bohb.minimizer[i]], marker = :star5, markersize = 8, color = :red)
                push!(plots, xy)
            end
        end
    end
    grid_plot = plot(plots..., layout=@layout([° _ _ _; ° ° _ _; ° ° ° _; ° ° ° °]))
    
    return grid_plot
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