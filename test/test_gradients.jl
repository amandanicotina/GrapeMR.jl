using Test
using GrapeMR

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = 0.0 
T1 = [1e8] 
T2 = [1e8] 
label = ["s1"]
target = ["min"]
spins = Spin(M0, T1, T2, offsets, ΔB1, target, label)

# Grape Parameters 
grape_params = GrapeParams(5, GrapeMR.spin_target, Dict("B1x" => true, "B1y" => true, "Bz" => false))

# Optimization Parameters
Tc, poly_start, poly_degree, max_iter = 0.5, 0.1, 1, 5;
opt_params = OptimizationParams(poly_start, poly_degree, max_iter);

# Parameters 
params = Parameters(grape_params, opt_params);

# Initial RF Pulse
B1ref = 1.0;
control_field = spline_RF(grape_params.N, Tc, B1ref)

# Dynamics
iso = dynamics(control_field, spins[1])
mag = iso.magnetization.dynamics
cost_vars = grape_params.cost_function(iso)
# Cost Variables
cost, adj_ini = cost_vars
# Adjoint Propagation
adj = backward_propagation(adj_ini, control_field, iso)
# True Gradient
true_grad_Bx = gradient(adj, mag, Ix)
true_grad_By = gradient(adj, mag, Iy)

# Finite difference
Δcf = 1e-6
fd_cf_Bx = finite_difference_field(spins[1], control_field, grape_params, "B1x", Δcf)
fd_cf_By = finite_difference_field(spins[1], control_field, grape_params, "B1y", Δcf)

using Plots
plot(fd_cf_Bx')
scatter!(true_grad_Bx')


# Test highest difference
tol = 1e-3
max_diff_Bx = [abs(true_grad_Bx[i] - fd_cf_Bx[i])/maximum([1, abs(true_grad_Bx[i]), abs(fd_cf_Bx[i])]) for i in eachindex(true_grad_Bx)]
max_diff_By = [abs(true_grad_By[i] - fd_cf_By[i])/maximum([1, abs(true_grad_By[i]), abs(fd_cf_By[i])]) for i in eachindex(true_grad_By)]

@test all(max_diff_Bx .< tol)
@test all(max_diff_By .< tol)









