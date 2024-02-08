using GrapeMR
using Plots

#spin = Spins([0.0; 0.0; 1.0], 0.5, 0.1, 0.0, "max")
spin = Spins([0.0; 0.0; 1.0], 1e8, 1e8, 0.0, "max")

N   = 500;
α   = π/2;
t_c = 1.0;
B1  = α/(γ_¹H*t_c); # [2πT]
B1x_arr = B1*ones(1, N); 
B1y_arr = zeros(1, N);  

init_control_field = InitialControlFields(N, B1x_arr, 0.0, B1y_arr, 0.0, t_c, 0.0, 0.0)

##### magnetization #####
mag = magnetization_ODE(init_control_field, spin)
plot(plot(mag[2:end, :]'))

##### OPTIMIZE #####
(iso, fields, cost_func) = grape_optimize(init_control_field, spin; max_iter=500, ϵ=1e-6)

##### PLOTS #####
data = [fields[1, :, i] for i in 1:101]
p = plot(data)
display(p)

