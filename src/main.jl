using GrapeMR
using ParameterSchedulers

##### INITIALIZATION #####
# Spins
M0 = [0.0; 0.0; 1.0];
T1 = [0.900, 0.070];
T2 = [0.300, 0.030];
target = ["min", "max"];
label = ["T1-900ms", "T1-70ms"];
# T1 = [1.830, 0.622, 2.430];
# T2 = [0.184, 0.092, 0.132];
# target = ["max", "min", "min"];# label = ["GD", "Yolk", "White"];


t_c = 0.5;
N   = 600;
B0  = [50.0];
# B0  = [-25.0, 0.0, 25.0];
# B0  = [-100.0, -75.0, -50.0, -25.0, 0.0, 25.0, 50.0, 75.0, 100.0]; # [Hz]
Bz  = zeros(1,N);
B1x = initial_field_spline(N, t_c)'; 
B1y = zeros(1,N);#initial_field_spline(N, t_c)';
ΔB1 = [1.0];

##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, B0, target, label, t_c, B1x, B1y, ΔB1, Bz);

##### OPTIMIZE #####
max_iter = 1000
lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, cost_saturation_contrast, [true false false]);
grape_output = @time grape(opt_params, field_init, spins, lr_scheduler; max_iter = max_iter,ϵ=1e-2); 

##### PLOTS #####
plot_magnetization_Mz_Mxy(grape_output.isochromats)
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

# TODO Fix lr_scheduler
# TODO Check if lr_scheduler is actually working and ϵ is not just constant
# TODO Improve how to choose cost function
# TODO Save Data


