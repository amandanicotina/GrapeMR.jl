using Test
using GrapeMR

# Parameters
N   = 2000;
t_c = 0.5; #[s]
cost_func = :saturation_contrast
grape_params = GrapeParams(N, cost_func, [true true false])

# Spin
M0  = [0.0, 0.0, 1.0];
T1  = [0.8, 0.4];
T2  = [0.6, 0.3];
B0  = [0.0];
ΔB1 = [1.0];
target = ["max", "min"];
label  = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

# Spin and RF objects
control_field = spline_RF(grape_params.N, t_c, 1.0)
spins         = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

# Dynamics
iso_vals = dynamics.(Ref(control_field), spins)

# Cost gradient
cost_vars = GrapeMR.cost_function.(iso_vals, grape_params.cost_function)
cost_grad = getindex.(cost_vars, 2)
adjoint   = backward_propagation.(control_field, iso_vals, cost_grad)

# True Gradient
t = range(0.0, control_field.t_control, length = grape_params.N)
Δt = round(t[2]-t[1], digits =5)
true_grad_Bx = GrapeMR.gradient.(adjoint, getfield.(getfield.(iso_vals, :magnetization), :dynamics), Ref(Ix), Δt)
true_grad_By = GrapeMR.gradient.(adjoint, getfield.(getfield.(iso_vals, :magnetization), :dynamics), Ref(Iy), Δt)

# Finite difference
Δcf = 1e-6 
fd_cf_Bx = finite_difference_field.(cost_func, control_field, spins, Δcf, "B1x")
fd_cf_By = finite_difference_field.(cost_func, control_field, spins, Δcf, "B1y")


# Test highest difference
tol = 1e-3
max_diff_Bx = maximum([maximum(abs.(true_grad_Bx[i] .- fd_cf_Bx[i])) for i in eachindex(true_grad_Bx)])
max_diff_By = maximum([maximum(abs.(true_grad_By[i] .- fd_cf_By[i])) for i in eachindex(true_grad_By)])

@test max_diff_Bx < tol 
@test max_diff_By < tol









