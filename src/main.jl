using GrapeMR

# Spin Parameters
<<<<<<< Updated upstream
M0      = [0.0, 0.0, 1.0]
ΔB1     = [1.0]
# offsets = -15:1:15
offsets = 1:1:1
T1 = [0.6, 0.1]
T2 = [0.3, 0.05]
label  = ["s1", "s2"]
target = ["max", "min"]
spins  = Spin(M0, T1, T2, offsets, ΔB1, target, label)

# Optimization Parameters
Tc = 1.0
poly_start  = 0.5
poly_degree = 1
max_iter    = 6000
opt_params  = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.euclidean_norm, Dict("B1x" => true, "B1y" => true, "Bz" => false))
=======
M0 = [0.0, 0.0, 1.0]
T1 = [0.6]
T2 = [0.3]
ΔB1 = [1.0]
B0 = 15.0
label  = ["s1"]
target = ["saturation"]
offset = -B0:1:B0
spins = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
max_iter = 6000
poly_start = 0.5
poly_degree = 1
opt_params = OptimizationParams(poly_start, poly_degree, max_iter)

# Grape Parameters 
N = 2000
cost = GrapeMR.euclidean_norm
fields2opt = Dict("B1x" => true, "B1y" => true, "Bz" => false)
grape_params = GrapeParams(N, cost, fields2opt)
>>>>>>> Stashed changes

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
<<<<<<< Updated upstream
B1ref = 1.0
control_field = hard_RF(grape_params.N, Tc, B1ref) 
=======
B1ref = 5.0
control_time = 0.5
control_field = spline_RF(N, control_time, B1ref)
>>>>>>> Stashed changes

# Run Optimization
grape_output = @time GrapeMR.grape(params, control_field, spins)

plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)
<<<<<<< Updated upstream
=======
plot_magnetization_time(grape_output.isochromats[2], grape_output.control_field.t_control)

# Cost Function Analysis
GrapeMR.run_cost_analysis(grape_output.control_field, spins[2], 50.0, 20, grape_params.cost_function)

# Save Optimization Data
folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/SimulationResults/"
save_grape_data(grape_output; folder_path)
>>>>>>> Stashed changes
