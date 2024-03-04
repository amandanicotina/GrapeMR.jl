using GrapeMR
using Plots

##### INITIALIZATION #####
# Spins
M0 = [0.0; 0.0; 1.0];
T1 = [0.100, 0.500];
T2 = [0.050, 0.300];
target = ["min", "max"];
label = ["-", "-"];
# T1 = [1.830, 0.622, 2.430];
# T2 = [0.184, 0.092, 0.132];
# target = ["max", "min", "min"];
# label = ["GD", "Yolk", "White"];

# Initial RF field
N   = 600;
αx  = π/2;
αy  = π/6;
t_c = 0.3;

time = range(0.0, t_c, N);
t    = time .- t_c/2;
rotx = rad2deg(αx)/360;
roty = rad2deg(αy)/360;
flip_x = rotx/diff(t)[1];
flip_y = roty/diff(t)[1];
BW_Hz = 300.0;
x     = BW_Hz.*t;
y     = BW_Hz.*t;
B0    = [-150.0, 0.0, 150.0] #Hz
B1x   = ((flip_x.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
B1y   = ((flip_y.*sinc.(y))./2π)';

##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, target, label, t_c, B1x, B1y, B0);

##### OPTIMIZE #####
opt_params = OptimizationParams(N, cost_saturation_contrast, [true true false]);
grape_output = @time grape(opt_params, field_init, spins; max_iter=25000, ϵ=1e-4); 

##### PLOTS #####
cost1 = grape_output.cost_values[1,:];
cost2 = grape_output.cost_values[2,:];
# cost3 = grape_output.cost_values[3,:];
# cost  = (cost1 + cost2 + cost3)
cost  = (cost1 + cost2)
cost[end] 
# cost = grape_output.cost_values[1,:];
plot_magnetization_Mz_Mxy(grape_output.isochromats[3:end])
plot_control_fields(grape_output.control_field) 
plot_control_fields(field_init) 
plot_magnetization_time(grape_output.isochromats[3], field_init.t_control)
plot_magnetization_time(grape_output.isochromats[4], field_init.t_control)
plot_cost_values(cost, opt_params)