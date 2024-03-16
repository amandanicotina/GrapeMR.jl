using GrapeMR
using Plots

##### INITIALIZATION #####
# Spins
M0 = [0.0; 0.0; 1.0];
# T1 = [0.180, 0.120];
# T2 = [0.120, 0.080];
T1 = [0.100, 0.500];
T2 = [0.050, 0.300];
target = ["max", "min"];
label = ["T1-100ms", "T2-500ms"];
# T1 = [1.830, 0.622, 2.430];
# T2 = [0.184, 0.092, 0.132];
# target = ["max", "min", "min"];# label = ["GD", "Yolk", "White"];

t_c = 0.6;
N   = 600;
#B0  = [0.0];
B0  = [-150.0, -100.0, -50, 0.0, 50.0, 100.0, 150.0]; # [Hz]
Bz  = zeros(1,N);
B1 = 1.0; # [Hz]
B1x = B1*initial_field_spline(N, t_c, B1)'; 
B1y = B1*initial_field_spline(N, t_c, B1)';
ΔB1 = [0.9, 1.0]

##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, B0, target, label, t_c, B1x, B1y, Bz);
plot_control_fields(field_init)

##### OPTIMIZE #####
opt_params = OptimizationParams(N, cost_saturation_contrast, [true true false]);
max_iter = 20000
# lr_scheduler = Step(start = 1e-1, decay=0.98, step_sizes=[(max_iter // 6) * 3, (max_iter // 6) * 2, max_iter // 6])
# lr_scheduler = Exp(start=1e-1, decay=0.99)
lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1)
grape_output = @time grape(opt_params, field_init, spins, lr_scheduler; max_iter=max_iter, ϵ=1e-4); 

##### PLOTS #####
plot_magnetization_Mz_Mxy(grape_output.isochromats)
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

##### Weigths & Biases #####
