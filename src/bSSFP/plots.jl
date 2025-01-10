theme(:default, 
    palette=:viridis,          # For discrete colors (lines, points)
    colorgradient=:viridis     # For continuous colors (colorbars, heatmaps)
)

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

"""
    plot_ss_offset_profile(ss_spin::Vector{GrapeMR.SteadyState})

Plot the bSSFP offset frequency profile showing magnitude and phase.

# Arguments
- `ss_spin`: Vector of steady-state spins at different offset frequencies

# Returns
- Tuple of three plots: (magnitude profile, phase profile, combined plot)

# Example
plots = plot_ss_offset_profile(ss_spins)
display(plots[1])  # Show
"""
function plot_ss_offset_profile(ss_spin::Vector{GrapeMR.SteadyState})
    colors = color_palette(10)
    data = calculate_steady_state(ss_spin)
    α = ss_spin[1].α

    # Plot 1: Magnitude Profile
    p1 = initialize_plot("bSSFP Offset Profile", "Offset Frequency (Hz)", "Magnitude")
        plot!(p1, data.offset, abs.(data.ss_matrix); label = "Matrix", color = colors[2], lw = 3, legend = :bottom)
        # plot!(p1, data.offset, abs.(data.ss_manual); label = "Manual", marker = :diamond, markersize=6, color = colors[4], lw = 3)
        scatter!(p1, data.offset, abs.(data.ss_geometric); label = "Geometric", marker = :star, markersize=6, color = colors[8], lw = 3)

    # Plot 2: Phase Profile
    p2 = initialize_plot("bSSFP Offset", "Offset Frequency (Hz)","Phase (rad)")
        yticks!([-π, 0, π], ["-π", "0", "π"])
        plot!(p2, data.offset, angle.(data.ss_matrix); label = "Matrix", color = colors[2], lw = 3)
        # plot!(p2, data.offset, angle.(data.ss_manual); label = "Manual", marker = :diamond, markersize=6, color = colors[4])
        scatter!(p2, data.offset, angle.(data.ss_geometric); label = "Geometric", marker = :star, markersize=6, color = colors[8], lw = 3)

    # Plot 4: Magnitude and Phase Profile
    α_deg = Int(round(rad2deg(α)))
    p3m = initialize_plot("bSSFP Offset Profile - α = $(α_deg)°", "", "Magnitude")
          plot!(p3m, data.offset, abs.(data.ss_matrix); label = false, color = colors[1], lw = 3.5)
    p3p = initialize_plot("", "Offset Frequency (Hz)", "Phase (rad)")
        yticks!([-π, -π/2, 0, π/2, π], ["-π", "-π/2", "0", "-π/2", "π"])
          plot!(p3p, data.offset, angle.(data.ss_matrix); label = false, color = colors[6], lw = 3.5)
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
    p4 = initialize_plot("bSSFP Offset Profile", "Offset Frequency (Hz)", "Magnitude")

    # Iterate over flip angles and add to plot
    colors = color_palette(flip_angles)
    for (i, (α, ss_sig)) in enumerate(flip_angle_profiles)
        α_deg = Int(round(rad2deg(α)))
        plot!(p4, offsets, abs.(ss_sig), label = "α = $(α_deg)°", lw = 3, color = colors[i], legend = :bottomleft)
    end

    return p4
end




function plot_magnetization_trajectory(ss_spin::Vector{GrapeMR.SteadyState})
    # Extract parameters
    offsets = collect(range(ss_spin[1].B0inho, ss_spin[end].B0inho, length(ss_spin)))
    
    # Extract Mxy and Mz components
    Mxy = steady_state_geometric.(ss_spin)  # Transverse magnetization
    Mz = steady_state_geometric_Mz.(ss_spin)  # Longitudinal magnetization
    
    # Create plot
    p = initialize_plot("bSSFP Magnetization", "Mxy", "Mz")
    
    # Plot geometric solution curve
    plot!(p, Mxy, Mz, 
          color=:black, 
          lw=1.5, 
          label="Ellipsoid")
    
    # Plot individual quivers with colors from colormap
    for i in 1:length(offsets)
        quiver!(p, [0], [0], 
                quiver=([Mxy[i]], [Mz[i]]), 
                lw=2.5, 
                label=false,          
                zcolor=offsets,
                colorbar_title="Offsets (Hz)",
                color=get(cgrad(:viridis, length(offsets)), i/length(offsets)))
    end
    
    p = add_arc_to_trajectory(p, Mxy, Mz, 1, 40)
    
    # Set plot limits
    xlims!(-1, 1)
    ylims!(0, 1)
    
    return p
end



function plot_bssfp_magnetization(Mxy, Mz, offsets, α; 
                                 show_arc::Bool=true, 
                                 arc_scale::Float64=0.5)
    # Create main plot
    p = plot(xlabel="Mxy", ylabel="Mz", 
            title="Magnetization - Offset Frequencies", 
            titlefontsize=15)
    
    # Plot geometric solution curve
    scatter!(p, Mxy, Mz, 
            color=:black, markersize=2, lw=0.6, label=false)
    
    # Add quivers with colors corresponding to the frequency offsets
    quiver!(p, zeros(length(offsets)), zeros(length(offsets)), 
            quiver=(Mxy, Mz), 
            lw=2.5, 
            label=false,
            zcolor=offsets,          
            colorbar_title="Frequency Offset (Hz)",
            color=:viridis)
    
    # Add flip angle arc if requested
    if show_arc
        # Create arc points
        function create_arc_points(start_vec, end_vec, α, n_points=50)
            θ1 = atan(start_vec[2], start_vec[1])
            θ2 = θ1 + α
            θ = range(θ1, θ2, length=n_points)
            radius = sqrt(start_vec[1]^2 + start_vec[2]^2)
            x = radius .* cos.(θ)
            y = radius .* sin.(θ)
            return x, y
        end
        
        # Choose vectors and create arc
        i = 1
        start_vec = [Mxy[i] * arc_scale, Mz[i] * arc_scale]
        end_vec = [Mxy[i+1] * arc_scale, Mz[i+1] * arc_scale]
        
        arc_x, arc_y = create_arc_points(start_vec, end_vec, α)
        plot!(arc_x, arc_y, 
              color=:black, 
              linestyle=:dash, 
              label=false)
        
        # Add angle label
        mid_point = floor(Int, length(arc_x)/2)
        annotate!(0.0, 0.5, 
                 text("α = $(round(rad2deg(α)))°", 12, :black))
    end
    
    # Set plot limits
    xlims!(-0.01, 0.5)
    ylims!(0, 0.7)
    
    return p
end

function plot_transverse_evolution(iso, spin; 
                                 time_unit=:ms,  # Options: :s, :ms, :μs
                                 show_ss=true)    # Option to show steady state line
    # Extract magnetization components
    Mx = iso.magnetization.dynamics[2, :]
    My = iso.magnetization.dynamics[3, :]
    Mxy = sqrt.(Mx.^2 + My.^2)
    
    # Create time array
    N = length(Mxy)
    t = range(0, iso.control_field.t_control, length=N)
    
    # Convert time based on chosen unit
    t_scaled = if time_unit == :ms
        t .* 1e3
        elseif time_unit == :μs
            t .* 1e6
        else
            t
    end
    
    # Get steady state value
    ss_value = abs(steady_state_geometric(spin))
    
    # Create plot
    p = plot(t_scaled, Mxy,
            xlabel="Time ($(time_unit))",
            ylabel="Transverse Magnetization |Mxy|",
            title="bSSFP Magnetization Evolution",
            label="Mxy",
            linewidth=2)
    
    # Add steady state line if requested
    if show_ss
        hline!([ss_value], 
               linestyle=:dash, 
               color=:red,
               label="Steady State",
               linewidth=1.5)
        
        # Add annotation for steady state value
        annotate!(maximum(t_scaled)-0.1*maximum(t_scaled), 
                 ss_value+0.02, 
                 text("SS = $(round(ss_value, digits=3))", 10, :red))
    end
    
    return p
end

"""
    add_arc_to_trajectory(p, Mxy, Mz, i1::Int, i2::Int)

Add an arc between two specified quivers on an existing magnetization trajectory plot.

# Arguments
- `p`: Existing plot
- `Mxy`: Vector of transverse magnetization values
- `Mz`: Vector of longitudinal magnetization values
- `i1`: Index of first vector
- `i2`: Index of second vector

# Returns
- Updated plot with arc
"""
function add_arc_to_trajectory(p, Mxy, Mz, i1::Int, i2::Int)
    # Calculate angles for the arc
    θ1 = atan(Mz[i1], Mxy[i1])
    θ2 = atan(Mz[i2], Mxy[i2])
    
    # Ensure we take the shorter path
    if abs(θ2 - θ1) > π
        if θ2 > θ1
            θ2 -= 2π
        else
            θ1 -= 2π
        end
    end
    
    # Create arc points
    n_points = 50
    θ = range(θ1, θ2, length=n_points)
    r = sqrt(Mxy[i1]^2 + Mz[i1]^2)  # Use radius of first vector
    
    arc_x = r .* cos.(θ)
    arc_y = r .* sin.(θ)
    
    # Add arc to plot
    plot!(p, arc_x, arc_y, 
          color=:black, 
          linestyle=:dash, 
          label=false,
          lw=1.5)
    
    # Add angle label at midpoint
    θ_mid = (θ1 + θ2)/2
    angle_deg = abs(rad2deg(θ2 - θ1))
    annotate!(p, 1.1*r*cos(θ_mid), 1.1*r*sin(θ_mid), 
             text("α = $(round(Int, angle_deg))°", 12, :black))
    
    return p
end

# Usage example:
# p = plot_magnetization_trajectory(ss_spin)
# Mxy = steady_state_geometric.(ss_spin)
# Mz = steady_state_geometric_Mz.(ss_spin)
# p = add_arc_to_trajectory(p, Mxy, Mz, 25, 75)  # Add arc between vectors 25 and 75


