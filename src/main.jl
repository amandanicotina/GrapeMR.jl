using GrapeMR
using Plots

spin = Spins([0.0; 1.0; 0.0], 0.5, 0.1, 0.0, "max")
#spin = Spins([0.5; 0.5; 0.0], 1e8, 1e8, 0.0, "max")

N   = 500;
α   = π/2;
t_c = 0.1;
B1  = α/(γ_¹H*t_c); # [2πT]
B1x_arr = B1*ones(1, N); 
B1y_arr = zeros(1, N);  

init_control_field = InitialControlFields(N, γ_¹H*B1x_arr, 0.0, B1y_arr, 0.0, t_c, 0.0, 0.0);

##### magnetization #####
mag = forward_propagation(init_control_field, spin);
iso = Magnetization((mag,),(spin,));
plot(plot(mag[2:end, :]'))
#mag_norm = forward_propagation(field_norm, spin_norm);
#iso_norm = Magnetization((mag_norm,),(spin_norm,));


##### OPTIMIZE #####
(iso_opt, fields_opt, cost_func) = grape_optimize(init_control_field, spin; max_iter=2000, ϵ=1e-4)

##### PLOTS #####
data = [fields_opt[1, :, i] for i in 1:101]
p = plot(data)
display(p)


p_mag = plot_magnetization(iso_opt, t_c)
display(p_mag)