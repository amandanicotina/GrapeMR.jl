
"""
    cost_offsets(control_field::ControlField, spin_system::Spin, offsets::Vector{Float64})

Calculates the cost values for a given control field (`control_field`) and spin system (`spin_system`) over a range of B0 inhomogeneities (`offsets`).

# Arguments
- `control_field::ControlField`: The control field object that contains the RF pulse sequence.
- `spin_system::Spin`: The spin system for which the cost function is to be calculated.
- `offsets::Vector{Float64}`: A vector of offset frequencies (in Hz) representing B0 inhomogeneities.

# Outputs
- `cost_values::Vector{Float64}`: A vector containing the cost values corresponding to each offset in `B0_values`.
"""
function cost_offsets(control_field::ControlField, spin_system::Spin, offsets::Vector{Float64}, cost::Symbol)
    ΔB1 = [1.0]
    spins_b0 = GrapeMR.Spin(spin_system.M_init, spin_system.T1, spin_system.T2, offsets, ΔB1, [spin_system.target], [spin_system.label])        
    isochromats = dynamics.(control_field, spins_b0)

    cost_terms = cost_function.(isochromats,  cost)  
    cost_profile = getindex.(cost_terms, 1) 

    return cost_profile
end

function cost_offsets(control_field::ControlField, spin_system::SteadyState, offsets::Vector{Float64})
    ΔB1 = [1.0]
    spins_b0 = GrapeMR.SteadyState(spin_system.M_init, spin_system.T1, spin_system.T2, offsets, ΔB1, [spin_system.target], 
                            [spin_system.label], spin_system.α, spin_system.Δϕ, spin_system.TR, spin_system.TE)      

    isochromats = dynamics.(control_field, spins_b0)

    cost_terms = cost_function.(isochromats, :spin_target)  
    cost_profile = getindex.(cost_terms, 1) 

    return cost_profile
end


"""
    create_cost_matrix(control_field::ControlField, spin_system::Spin, offsets::Vector{Float64}, b1_inhomogeneities::Vector{Float64})

Generates a cost matrix representing the cost function values across a range of B0 and B1 inhomogeneities.

# Arguments
- `control_field::ControlField`: The control field object that contains the RF pulse sequence
- `spin_system::Spin`: The spin system for which the cost function is to be calculated.
- `offsets::Vector{Float64}`: A vector of offset frequencies (in Hz) representing B0 inhomogeneities.
- `b1_inhomogeneities::Vector{Float64}`: A vector representing B1 inhomogeneity percentages.

# Outputs
- `cost_matrix::Matrix{Float64}`: A 2D matrix where each element represents the cost value for a particular 
                                        combination of B0 and B1 values. The matrix is normalized by its maximum value.
"""
function create_cost_matrix(control_field::ControlField, spin_system::Spin, offsets::Vector{Float64}, b1_inhomogeneities::Vector{Float64}, cost::Symbol)
    cost_matrix = Matrix{Float64}(undef, length(b1_inhomogeneities), length(offsets))

    for (i, b1) in enumerate(b1_inhomogeneities)
        spins_b0 = GrapeMR.Spin(spin_system.M_init, spin_system.T1, spin_system.T2, offsets, b1, [spin_system.target], [spin_system.label])        
        isochromats = dynamics.(control_field, spins_b0)

        cost_terms = cost_function.(isochromats, cost)  
        cost_profile = getindex.(cost_terms, 1)  
        cost_matrix[i,:] = cost_profile
    end

    cost_matrix = cost_matrix ./ maximum(cost_matrix)
    return cost_matrix
end
function create_cost_matrix(control_field::ControlField, spin_system::SteadyState, offsets::Vector{Float64}, b1_inhomogeneities::Vector{Float64}, cost::Symbol)
    cost_matrix = Matrix{Float64}(undef, length(b1_inhomogeneities), length(offsets))

    for (i, b1) in enumerate(b1_inhomogeneities)
        spins_b0 = GrapeMR.SteadyState(spin_system.M_init, spin_system.T1, spin_system.T2, offsets, b1, [spin_system.target], 
                                        [spin_system.label], spin_system.α, spin_system.Δϕ, spin_system.TR, spin_system.TE)     
        isochromats = dynamics.(control_field, spins_b0)

        cost_terms = cost_function.(isochromats, cost)  
        cost_profile = getindex.(cost_terms, 1)  
        cost_matrix[i,:] = cost_profile
    end

    cost_matrix = cost_matrix ./ maximum(cost_matrix)
    return cost_matrix
end


"""
    plot_cost_offset(cost_profile::Vector{Float64}, B0_values::Vector{Float64})

Creates a plot of the cost values as a function of B0 offset frequencies.

# Arguments
- `cost_profile::Vector{Float64}`: A vector containing the cost values to be plotted.
- `B0_values::Vector{Float64}`: A vector of offset frequencies (in Hz) corresponding to the cost values.

# Outputs
- `p::Plot`: A plot object displaying the cost function offset profile.
"""
function plot_cost_offset(cost_profile::Vector{Float64}, offsets::Vector{Float64})
    p = plot(xlabel = "Offset [Hz]",
             ylabel = "Cost Value",
             title  = "Cost Function Offset Profile",
             guidefontsize = 12, 
             legendfontsize = 10,
             tickfontsize = 10,
             titlefontsize = 12,
             framestyle = :box, 
             grid = false)
    plot!(p, offsets, cost_profile, label = "min", lw = 2)
    return p
end



"""
    heatmap_cost(cost_matrix::Matrix{Float64}, offsets::Vector{Float64}, b1_inhomogeneities::Vector{Float64})

Generates a heatmap plot of the cost function values over a range of B0 and B1 inhomogeneities.

# Arguments
- `cost_matrix::Matrix{Float64}`: A 2D matrix representing the cost values for combinations of B0 and B1 inhomogeneities.
- `offsets::Vector{Float64}`: A vector of offset frequencies (in Hz) used to label the x-axis.
- `b1_inhomogeneities::Vector{Float64}`: A vector representing B1 inhomogeneity percentages used to label the y-axis.

# Outputs
- `h::Plot`: A heatmap plot showing the cost function map with a color bar indicating cost values.
"""
function heatmap_cost(cost_matrix::Matrix{Float64}, offsets::Vector{Float64}, b1_inhomogeneities::Vector{Float64})
    heatmap(offsets, b1_inhomogeneities, cost_matrix, color=:viridis, framestyle=:box,
            guidefontsize = 12, 
            legendfontsize = 10,
            tickfontsize = 10,
            titlefontsize = 12,
            xlabel="Offset [Hz]", 
            ylabel="B1 Inhomogeneity [%]", 
            title="Cost Function Map",
            colorbar_title="Cost Value")
end


"""
    countour_cost(cost_matrix::Matrix{Float64})

Creates a contour plot of the cost function values over the range of B0 and B1 inhomogeneities.

# Arguments
- `cost_matrix::Matrix{Float64}`: A 2D matrix representing the cost values for combinations of B0 and B1 inhomogeneities.

# Outputs
- `c::Plot`: A contour plot of the cost function map.
"""
function countour_cost(cost_matrix::Matrix{Float64})
    contourf(cost_matrix, color=:viridis, framestyle=:box,
             guidefontsize = 12, 
             legendfontsize = 10,
             tickfontsize = 10,
             titlefontsize = 12,
             xlabel="Offset [Hz]", 
             ylabel="B1 Inhomogeneity [%]", 
             title="Cost Function Map", 
             colorbar_title="Cost Value",
             levels=20)
end



"""
    run_cost_analysis(grape_output::GrapeOutput, offset::Float64, b1_inhomogeneity_percent::Int)

Runs a complete cost analysis for the specified `GrapeOutput` object over a range of B0 and B1 inhomogeneities, 
and generates plots to visualize the cost function.

# Arguments
- `control_field::ControlField`: The control field object that contains the RF pulse sequence.
- `spin_system::Spin`: The spin system for which the cost function is to be calculated.
- `offset::Float64`: The range of offset frequencies (in Hz) for the B0 inhomogeneities.
- `b1_inhomogeneity_percent::Int`: The percentage of B1 inhomogeneity to consider for the analysis.

# Workflow
1. Generates a range of B0 values from `-offset` to `offset`.
2. Computes cost values for the specified B0 and B1 inhomogeneities using `cost_offsets` and `create_cost_matrix`.
3. Plots the cost function offset profile using `plot_cost_offset`.
4. Generates a heatmap of the cost matrix using `heatmap_cost`.
5. Creates a contour plot of the cost matrix using `countour_cost`.

# Example
run_cost_analysis(grape_output, 50.0, 30)
"""
function run_cost_analysis(control_field::ControlField, spin_system::Spins, offset::Float64, b1_inhomogeneity_percent::Int, cost::Symbol)
# It has to somehow be grape_output because I need the cost func information
    offsets = collect(-offset:1:offset)
    b1_variation = b1_inhomogeneity_percent / 100
    b1_range = collect(range(1.0 - b1_variation, 1.0 + b1_variation, length=5))

    cost_profile = cost_offsets(control_field, spin_system, offsets,  cost)
    cost_matrix = create_cost_matrix(control_field, spin_system, offsets, b1_range, cost)

    p1 = plot_cost_offset(cost_profile, offsets)
    h1 = heatmap_cost(cost_matrix, offsets, b1_range)
    c1 = countour_cost(cost_matrix)
    
    display(p1)
    display(h1)
    display(c1)
end

# # run_cost_analysis(grape_output, 50.0, 30)
