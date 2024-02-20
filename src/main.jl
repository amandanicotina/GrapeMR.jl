using GrapeMR
using Plots

# Create abstract arrays for T1 etc for N spins and use as input of the normalization function

##### SPINS #####
M0 = [1.0; 0.0; 0.0];
T1 = 0.5;
T2 = 0.25;

##### INITIAL RF FIELD #####
N   = 500;
α   = π/2;
t_c = 0.5;
B1  = α/t_c;
B1x, B1y = B1*ones(Float64, 1, N), B1*ones(Float64, 1, N)
B0, ΔB0  = zeros(1, N), zeros(1, N);
 
(spins, init_control_field) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

##### INITIAL MAGNETIZATION #####
mag  = forward_propagation(init_control_field, spins);
iso  = Magnetization((mag,),(spins,));
back = backward_propagation(init_control_field, iso, "Target One Spin")

##### OPTIMIZE #####
opt_params = OptimizationParams(N, "Target One Spin", [true true false]);
(iso_opt, fields_opt) = grape_optimize(opt_params, init_control_field, spins; max_iter=1000, ϵ=1e-4); 

##### PLOTS #####
p_mag = plot_magnetization(iso_opt, t_c)
p_ini = plot_magnetization(iso, t_c)
p_field = plot_control_fields(fields_opt)
p_init = plot_control_fields(init_control_field)
display(p_mag)
display(p_field)