# using Plots
# using Base.Threads

# # Function to initialize spins with offset and B1 values
# function initialize_spins(M0, B0_vals, B1_vals, targetoff, labeloff)
#     return GrapeMR.Spin(M0, [0.3], [0.08], B0_vals, B1_vals, targetoff, labeloff)
# end

# # Function to calculate magnetization dynamics and cost function for each spin
# function calculate_iso_and_cost(spins, control_field)
#     isos = Isochromat[]
#     costs = []

#     for spin in spins
#         mag = forward_propagation(control_field, spin)
#         dyn = GrapeMR.Magnetization(mag)
#         iso = Isochromat(dyn, spin)
#         push!(isos, iso)
#         push!(costs, cost_function(iso, :target_one_spin))
#     end

#     return isos, costs
# end

# # Function to create the cost matrix for B0 and B1 values
# function create_cost_matrix(B0_values, B1_values, control_field, M0, targetoff, labeloff)
#     cost_matrix = zeros(length(B1_values), length(B0_values))

#     for (i, b1) in enumerate(B1_values)
#         for (j, b0) in enumerate(B0_values)
#             spin = GrapeMR.Spin(M0, [0.3], [0.08], b0, b1, targetoff, labeloff)
#             mag = forward_propagation(control_field, spin[1])
#             dyn = GrapeMR.Magnetization(mag)
#             iso = Isochromat(dyn, spin[1])
#             cost_func = cost_function(iso, :target_one_spin)

#             cost_matrix[i, j] = cost_func[1][1]
#         end
#     end
#     return cost_matrix
# end

# # Function to plot the cost function heatmap
# function plot_heatmap(cost_matrix, B0_values, B1_values)
#     heatmap(cost_matrix,
#         xlabel="B0 Inhomogeneity",
#         ylabel="B1 Inhomogeneity",
#         title="Cost Function Map",
#         # xticks=(1:length(B0_values), round.(B0_values, digits=2)),
#         # yticks=(1:length(B1_values), round.(B1_values, digits=2))
#     )
# end

# # Main function to run the whole process
# function run_simulation()
#     # Parameters
#     labeloff  = ["T1 = 300ms"]
#     targetoff = ["min"]
#     offset = 10  # [Hz]
#     offset_vals = collect(-offset:2:offset)
#     n = 5
#     B1_vals = collect(range(0.7, stop=1.3, length=n))

#     # Initialize spins
#     spins_offset = initialize_spins(M0, offset_vals, [1.0], targetoff, labeloff)
#     spins_offset_B1 = initialize_spins(M0, offset_vals, B1_vals, targetoff, labeloff)

#     # Calculate magnetization dynamics and cost function for each spin
#     isos, costs = calculate_iso_and_cost(spins_offset, grape_output.control_field)

#     # Plot offset cost (can be customized based on specific needs)
#     p1 = plot_cost_offset(isos, :target_one_spin)
#     display(p1)

#     # Extract B0 and B1 inhomogeneity values for heatmap
#     B0_values = [spin.B0inho for spin in spins_offset_B1]
#     B1_values = [spin.B1inho for spin in spins_offset_B1]

#     # Create cost matrix
#     @time cost_matrix = create_cost_matrix(B0_values, B1_values, grape_output.control_field, M0, targetoff, labeloff)

#     # Plot heatmap
#     h1 = plot_heatmap(cost_matrix, B0_values, B1_values)
#     display(h1)
# end

# # Example threading implementation for parallel computation
# function create_cost_matrix_parallel(B0_values, B1_values, control_field, M0, targetoff, labeloff)
#     cost_matrix = zeros(length(B1_values), length(B0_values))

#     @time @threads for i in enumerate(B1_values)
#         b1 = B1_values[i]

#         for j in enumerate(B0_values)
#             b0 = B0_values[j]

#             spin = GrapeMR.Spin(M0, [0.3], [0.08], b0, b1, targetoff, labeloff)

#             mag = forward_propagation(control_field, spin)
#             dyn = GrapeMR.Magnetization(mag)
#             iso = Isochromat(dyn, spin)
#             cost_func = cost_function(iso, :target_one_spin)

#             cost_matrix[i, j] = cost_func[1][1]
#         end
#     end

#     return cost_matrix
# end

# # Call the main function to run the simulation
# run_simulation()
