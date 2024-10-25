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

# Optimization Parameters
Tc = 0.5
poly_start = 0.1
poly_degree = 1.0
max_iter = 100.0
opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, "spin_target", [true true false])

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref) 

# Run Optimization
grape_output = @time grape(params, control_field, spins); 