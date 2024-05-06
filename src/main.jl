using GrapeMR
using ParameterSchedulers

##### INITIALIZATION #####
# Spins
M0 = [1.0; 0.0; 0.0];
# T1 = [0.900, 0.070];
# T2 = [0.300, 0.030];
# target = ["min", "max"];
# label = ["T1-900ms", "T1-70ms"];
# T1 = [1.830, 0.622, 2.430];
# T2 = [0.184, 0.092, 0.132];
# label = ["GD", "Yolk", "White"];
# target = ["max", "min", "min"];

T1 = [0.5];
T2 = [0.25];
target = ["max"];
label = ["T1-500ms"];

t_c = 0.2;
N   = 600;
B0  = [0.0];
Bz  = zeros(1,N);
B1x = 1*initial_field_spline(N, t_c)'; 
B1y = zeros(1,N);#initial_field_spline(N, t_c)';
ΔB1 = [1.0];

##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, B0, target, label, t_c, B1x, B1y, ΔB1, Bz);
plot_control_fields(field_init)

##### OPTIMIZE #####
max_iter = 5000
lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, cost_target_one_spin, [true false false]);
grape_output = @time grape(opt_params, field_init, spins, lr_scheduler; max_iter = max_iter,ϵ=1e-2); 

##### PLOTS #####
plot_magnetization_Mz_Mxy(grape_output.isochromats)
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

# TODO Fix lr_scheduler
# TODO Check if lr_scheduler is actually working and ϵ is not just constant
# TODO Save Data



