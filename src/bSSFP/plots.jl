function calculate_steady_state(ss_spin::Vector{GrapeMR.SteadyState})
    offset     = collect(range(ss_spin[1].B0inho, ss_spin[end].B0inho, length(ss_spin)))
    ss_mat_vec = steady_state_matrix.(ss_spin)
    ss_man_vec = steady_state.(ss_spin)
    ss_geo_vec = steady_state_geometric.(ss_spin)

    SteadyStateData(
        offset,
        ss_mat_vec,
        ss_man_vec,
        ss_geo_vec
    )
end

function calculate_steady_state(ss_spin::Vector{GrapeMR.SteadyState}, flip_angles::Vector{Float64})
    offset = collect(range(ss_spin[1].B0inho, ss_spin[end].B0inho, length(ss_spin)))
    flip_angle_profiles = Vector{Tuple{Float64, Vector{ComplexF64}}}()
    s = ss_spin[1];
    for flip_angle ∈ flip_angles
        ss_spin_flip = GrapeMR.SteadyState(s.M_init, s.T1, s.T2, offset, s.B1inho, [s.target], [s.label], flip_angle, s.Δϕ, s.TR, s.TE)
        ss_mat_flip  = steady_state_matrix.(ss_spin_flip)
        push!(flip_angle_profiles, (flip_angle, ss_mat_flip))
    end
    return offset, flip_angle_profiles
end

function plot_ss_offset_profile(ss_spin::Vector{GrapeMR.SteadyState})
    colors = color_palette(1:7)
    data = calculate_steady_state(ss_spin)

    # Plot 1: Magnitude Profile
    p1 = plot(xlabel = "Offset Frequency (Hz)", ylabel = "Magnitude", title = "bSSFP Off-Resonance Profile", 
             titlefontsize = 15, framestyle = :box, grid = false)
        plot!(p1, data.offset, abs.(data.ss_matrix); label = "Matrix", color = colors[1], lw = 2)
        plot!(p1, data.offset, abs.(data.ss_manual); label = "Manual", marker = :diamond, color = colors[4], lw = 2)
        scatter!(p1, data.offset, abs.(data.ss_geometric); label = "Geometric", marker = :star, color = colors[7], lw = 2)

    # Plot 2: Phase Profile
    p2 = plot(xlabel = "Offset Frequency (Hz)", ylabel = "Phase (rad)", yticks = ([-π, 0, π], ["-π", "0", "π"]),
              title = "bSSFP Off-Resonance Profile", titlefontsize = 15, framestyle = :box, grid = false)
        plot!(p2, data.offset, angle.(data.ss_matrix); label = "Matrix", color = colors[1], lw = 2)
        plot!(p2, data.offset, angle.(data.ss_manual); label = "Manual", marker = :diamond, color = colors[4])
        scatter!(p2, data.offset, angle.(data.ss_geometric); label = "Geometric", marker = :star, color = colors[7], lw = 2)

    # Plot 4: Magnitude and Phase Profile
    p3m = plot(ylabel = "Magnitude", title = "bSSFP Off-Resonance Profile", titlefontsize = 15, framestyle = :box, grid = false)
          plot!(p3m, data.offset, abs.(data.ss_matrix); label = "Matrix", color = colors[1], lw = 2)
    p3p = plot(xlabel = "Offset Frequency (Hz)", ylabel = "Phase (rad)", yticks = ([-π/2, 0, -π/2], ["-π/2", "0", "-π/2"]), framestyle = :box, grid = false)
          plot!(p3p, data.offset, angle.(data.ss_matrix); label = "Matrix", color = colors[6], lw = 2)
    p3  = plot(p3m, p3p, layout = (2, 1))

    # Display Plots
    display(p1)
    display(p2)
    display(p3)

    return (p1, p2, p3)
end

function plot_ss_flip_angle_profile(ss_spin::Vector{GrapeMR.SteadyState}, flip_angles::Vector{Float64})

    (offsets, flip_angle_profiles) = calculate_steady_state(ss_spin, flip_angles);
    # Plot flip angle profiles
    p4 = plot(xlabel = "Offset Frequency (Hz)", ylabel = "Magnitude", title = "bSSFP Off-Resonance Profile",
              guidefontsize = 12, legendfontsize = 10, tickfontsize = 10, titlefontsize = 12, framestyle = :box, grid = false)

    # Iterate over flip angles and add to plot
    colors = color_palette(flip_angles)
    for (i, (α, ss_sig)) in enumerate(flip_angle_profiles)
        α_deg = Int(round(rad2deg(α)))
        plot!(p4, offsets, abs.(ss_sig), label = "α = $(α_deg)°", lw = 2, color = colors[i], legend = :bottomleft)
    end

    return p4
end


# function long_mag()
#     # Plot 3: Longitudinal Magnetization
#     p3 = plot(ss_data.offset, abs.(ss_data.ss_matrix[:,3]); label = "Matrix", marker = :diamond, color = colors[3], lw = 2,
#     xlabel = "Resonant Frequency (Hz)", ylabel = "Mz - Magnetization", title = "bSSFP Off-Resonance Profile", titlefontsize = 15)
#     scatter!(p3, ss_data.offset, ss_data.ss_geometric[:,2]; label = "Geometric", marker = :circle, color = colors[6], lw = 2)
# end
