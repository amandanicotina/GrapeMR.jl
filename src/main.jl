using GrapeMR

# Spin System
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = 0.0
offsets = collect(-B0:1:B0)

# Water
T1_water = 0.5
T2_water = 0.1
label_water = "S1"
target_water = "[0.0, 1.0, 0.0]"

spins = GrapeMR.Spin(M0, [T1_water], [T2_water], offsets, ΔB1, [target_water], [label_water])
# Grape Parameters 
grape_params = GrapeParams(1500, "spin_target", [true true false])

# Optimization Parameters
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.01, 1.0, 15), 1500)
# plot_bohb(bohb)

# spin = grape_output.isochromats[1].spin
# @time run_cost_analysis(grape_output.control_field, spin, 50.0, 50, grape_params.cost_function)

# Save data
# folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
# experiment_folder = save_grape_data(grape_output; folder_path)
# go = load_grape_data(experiment_folder)

# Export data
# export_bruker(grape_output)

# # Plots
# plot_cost_values(threads_grape_output.cost_values, grape_params)
# plot_magnetization_2D(grape_output.isochromats)
# plot_magnetization_control_field(threads_grape_output.control_field, grape_output.isochromats)
# plot_control_fields(grape_output.control_field; unit="Hz")
# plot_magnetization_time(grape_output.isochromats[11], grape_output.control_field.t_control)