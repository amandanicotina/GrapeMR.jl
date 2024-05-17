using GrapeMR
using ParameterSchedulers

##### INITIALIZATION #####
# Parameters
N   = 500;
t_c = 0.5; #[s]

# Spin
M0  = [0.0; 0.0; 1.0];
T1  = [0.5]*1e8;
T2  = [0.3]*1e8;
B0  = [-50.0, -25.0, 0.0, 25.0, 50.0];
# B0  = [0.0];
target = ["max"];
label  = ["T1-$(Int(T1[1]*1e3))ms"];
# T1 = [1.830, 0.622, 2.430];
# T2 = [0.184, 0.092, 0.132];
# label = ["GD", "Yolk", "White"];
# target = ["max", "min", "min"];

# RFs
ΔB1, Bz = [1.1, 1.0, 0.9], zeros(1,N), 0.0;
B1x = initial_field_spline(N, t_c)'; # (1, N); # initial_field_spline(N, t_c)'; # 
B1y = ones(1,N); # initial_field_spline(N, t_c)'; # 

# Spin and RF objects
control_field = ControlField(B1x, B1y, 1.0, Bz, t_c)
spins  = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)
plot_control_fields(control_field)

##### OPTIMIZE #####
max_iter     = 800
lr_scheduler = Poly(start=1e-2, degree=3, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, :target_one_spin, [true true false]);
grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 

##### PLOTS #####
plot_magnetization_Mz_Mxy(grape_output.isochromats)
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

##### SAVE DATA #####
folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
GrapeMR.save_grape_data(grape_output; folder_path)


# TODO Fix lr_scheduler
# TODO Check if lr_scheduler is actually working and ϵ is not just constant
# TODO Save Data

