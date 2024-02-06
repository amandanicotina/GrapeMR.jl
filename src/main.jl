using GrapeMR
using Plots

spin1 = Spins([0.5; 0.5; 0.0], 0.5, 0.1, 0.0, "max")

N   = 500;
α   = π/2;
t_c = 1.0;
B1  = α/(γ_¹H*t_c); # [2πT]
B1x_arr = B1*ones(1, N); 
B1y_arr = zeros(1, N);  

init_control_field = InitialControlFields(N, B1x_arr, 0.0, B1y_arr, 0.0, t_c, 0.0, 0.0)

##### OPTIMIZE #####
(iso, fields, cost_func) = grape_optimize(init_control_field, spin1; max_iter=500, ϵ=1e-6)

##### PLOTS #####
data = [fields[1, :, i] for i in 1:101]
p = plot(data)
display(p)

ploting = false
if ploting
    mag_ini = magnetization_ODE(init_control_field, spin1)
    isoInit = Magnetization((mag_ini,), (spin1,))
    (pi, p2) = plot_magnetization(isoInit, 10.0, t_c)
    (display(pi))   
    (p1, p2) = plot_magnetization(iso, 10.0)
    display(p1)
    mag_1 = magnetization_ODE(init_control_field, spin1)
    time  = range(0.0, t_c, length=N+1)
    plot(time, mag_1[2:end,:]') 
end
