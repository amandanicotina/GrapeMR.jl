using GrapeMR

# Spin Parameters
M0     = [0.0, 0.0, 1.0]
ΔB1    = [1.0]
# offset = -1:0.1:1
offset = 1:1:1
T1 = [1.35]
T2 = [0.05]
label  = ["Blood"]
target = ["-"]
spins  = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters
Tc = 0.1
poly_start  = 0.1
poly_degree = 1.0
max_iter    = 1
opt_params  = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.euclidean_norm, [true true false])

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref) 

# Run Optimization
# @time GrapeMR.grape(params, control_field, spins); 
# @profview GrapeMR.grape(params, control_field, spins); 

# Metal tests
out = grape(params, control_field, spins);
grape_output = @time metal_grape(params, control_field, spins);
# grape_output = @time threads_metal_grape(params, control_field, spins); 
