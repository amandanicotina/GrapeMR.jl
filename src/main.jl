using GrapeMR
using Plots

##### SPINS #####
spin = Spins([0.0; 0.0; 1.0], 0.5, 0.1, 0.0, "max")
#spin = Spins([0.5; 0.5; 0.0], 1e8, 1e8, 0.0, "max")

##### INITIAL RF FIELD #####
N    = 500;
α    = π/2;
t_c  = 0.05;
B1x_max, B1y_max = 2π/t_c, 2π/t_c;

time = range(0.0, t_c, N) 
t = time .- t_c/2;
flip_rads = α/diff(t)[1]
BW_Hz = 500;
x = 2*BW_Hz.*t/2
y = 2*BW_Hz.*t/2
B1x_arr = flip_rads*sinc.(x); # sinc(x) = sin(πx)/(πx)
B1y_arr = flip_rads*sinc.(y); #;  

B0  = zeros(1, N);
ΔB0 = zeros(1, N);

init_control_field = ControlFields((B1x_arr./B1x_max)', (B1y_arr./B1x_max)', 
                                    B1x_max, B1y_max, t_c, B0, ΔB0); 

##### INITIAL MAGNETIZATION #####
mag = forward_propagation(init_control_field, spin);
iso = Magnetization((mag,),(spin,));
back = backward_propagation(init_control_field, iso, "Target One Spin")

g = get_gradient(init_control_field, iso, Iy, "Target One Spin")
b = update_control_field(init_control_field, iso, Ix, "Target One Spin", 1e-6)
plot(b')

##### OPTIMIZE #####
opt_params = OptimizationParams(N, "Target One Spin", [true true false]);
(iso_opt, fields_opt) = grape_optimize(opt_params, init_control_field,  spin; max_iter=2000, ϵ=1e-5); 

##### PLOTS #####
p_mag = plot_magnetization(iso_opt, t_c)
p_ini = plot_magnetization(iso, t_c)
p_field = plot_control_fields(fields_opt)
p_init = plot_control_fields(init_control_field)
display(p_mag)
display(p_field)