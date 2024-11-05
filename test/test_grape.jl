using GrapeMR, Test

# Tolerance to past the tests
tol = 1e-2

#########################################################
#                Cost: Euclidean Norm                   #
#########################################################
# General Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = [0.0]

# Spins Parameters
T1 = [0.3]
T2 = [0.08]
label = ["S1"]
target = ["-"]

spins = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.euclidean_norm, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
Tc = 1.0
poly_start = 0.1
poly_degree = 1
max_iter = 1000
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output_en = @time grape(params, control_field, spins);

@test isapprox(round(grape_output_en.cost_values[end], digits=2), 0; atol=tol)

#########################################################
#              Cost: Saturation Contrast                #
#########################################################
# General Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = [0.0]

# Spins Parameters
T1 = [1.0, 0.1]
T2 = [0.6, 0.05]
label = ["S1", "S2"]
target = ["min", "max"]
spins = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.saturation_contrast, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
Tc = 1.0
poly_start = 0.1
poly_degree = 1
max_iter = 500
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output_sc = @time grape(params, control_field, spins);

@test isapprox(round(grape_output_sc.cost_values[end], digits=2), 0; atol=tol)


#########################################################
#    Cost: Euclidean Norm    CF: Sinc Function          #
#########################################################
# General Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = [0.0]

# Spins Parameters
T1 = [0.3]
T2 = [0.08]
label = ["S1"]
target = ["-"]

spins = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.euclidean_norm, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
Tc = 1.0
poly_start = 0.5
poly_degree = 1
max_iter = 500
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = sinc_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output = @time grape(params, control_field, spins);

@test isapprox(round(grape_output.cost_values[end], digits=2), 0; atol=tol)


#########################################################
#    Cost: Euclidean Norm    CF: Hard Pulse             #
#########################################################
# General Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = [0.0]

# Spins Parameters
T1 = [0.3]
T2 = [0.08]
label = ["S1"]
target = ["-"]

spins = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.euclidean_norm, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
Tc = 1.0
poly_start = 0.5
poly_degree = 1
max_iter = 250
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = hard_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output = @time grape(params, control_field, spins);

@test isapprox(round(grape_output.cost_values[end], digits=2), 0; atol=tol)



#########################################################
#    Cost: Euclidean Norm    CF: Gaussian              #
#########################################################
# General Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = [0.0]

# Spins Parameters
T1 = [1e8]
T2 = [1e8]
label = ["S1"]
target = ["-"]

spins = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(1500, GrapeMR.spin_target, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
Tc = 0.05
poly_start = 0.5
poly_degree = 2
max_iter = 100
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Parameters 
params = Parameters(grape_params, opt_params)

# Initial RF Pulse
B1ref = 1.0
control_field = gaussian_RF(grape_params.N, Tc, B1ref)

# Run Optimization
grape_output = @time grape(params, control_field, spins)

plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
plot_cost_values(grape_output.cost_values, grape_output.params.grape_params)
grape_output.cost_values[end]

@test isapprox(round(grape_output.cost_values[end], digits=2), 0; atol=tol)

