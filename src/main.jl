using GrapeMR


spin1 = Spins([1.0; 0.0; 0.0], 0.5, 0.3, 0.0, "max")
spin2 = Spins([0.5; -0.5; 0.0], 0.1, 0.05, 0.0, "min")

N   = 1000;
α   = π/2;
t_c = 2.5;
B1  = α/(γ_¹H*t_c);
B1_arr = B1*ones(1, N);  
#x = range(-4,4,length=N)
#B1_arr = sin.(α*x)./(α*x)

init_control_field = InitialControlFields(N, transpose(B1_arr), 1.0, 0.0, t_c, 0.0, 0.0)
plot(init_control_field.B1_initial_control)
##### OPTIMIZE #####
(iso, fields) = grape_optimize(init_control_field, spin1; ϵ=1e-5)



##### PLOTS #####
data = [fields[1, :, i] for i in 1:101]
p = plot(data)
display(p)

ploting = true
if ploting
    mag_ini = magnetization_ODE(init_control_field, spin1)
    isoInit = Magnetization((mag_ini,), (spin1,))
    (pi, p2) = plot_magnetization(isoInit, 10.0, t_c)
    (display(pi))   
    (p1, p2) = plot_magnetization(iso, 10.0)
    display(p1)
    mag_1 = magnetization_ODE(init_control_field, spin1)
    time = range(0.0, t_c, length=N+1)
    plot(time, mag_1[2:end,:]') 
end
