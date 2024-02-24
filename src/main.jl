using GrapeMR
using Plots

# Create abstract arrays for T1 etc for N spins and use as input of the normalization function

##### SPINS #####
M0 = [0.0; 0.0; 1.0];
T1 = [0.5, 0.7];
T2 = [0.2, 0.3];

DEBUG = false

# Initial RF field
#N   = 500;
#t_c = 0.05;
#B1  = 10.0;
#B1x, B1y = B1*ones(Float64, 1, N), ones(Float64, 1, N)
#B0, ΔB0  = zeros(1, N), zeros(1, N);

N   = 500;
αx  = π/2;
αy  = π/6;
t_c = 0.1;

time = range(0.0, t_c, N);
t    = time .- t_c/2;
rotx = rad2deg(αx)/360;
roty = rad2deg(αy)/360;

flip_x = rotx/diff(t)[1];
flip_y = roty/diff(t)[1];

BW_Hz = 500.0;
x     = BW_Hz.*t;
y     = BW_Hz.*t;
B0    = zeros(1, N);
#B1x   = ((flip_x.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
#B1y   = zeros(1, N);

B1x   = ((flip_x.*sinc.(x))./2π)';
B1y   = ((flip_y.*sinc.(y))./2π)';

(spins, field_init) = normalization(M0, t_c, T1, T2, B1x, B1y, B0)

if DEBUG
    for spin ∈ spins
        ##### INITIAL MAGNETIZATION #####
        mag  = forward_propagation(field_init, spin);
        iso  = Magnetization((mag,),(spin,));
        back = backward_propagation(field_init, iso, "Target One Spin")
        plot(back')
        ##### GRADIENTS #####
        g = get_gradient(field_init, iso, Iy, "Target One Spin");
        b = update_control_field(field_init, iso, Ix, "Target One Spin", 1e-3)
        plot(g')
    end
end

##### OPTIMIZE #####
opt_params = OptimizationParams(N, "Target One Spin", [true true false]);
grape_output = grape_optimize(opt_params, field_init, spins; max_iter=1000, ϵ=1e-3); 

##### PLOTS #####
for (idx, (iso_opt, fields_opt)) ∈ enumerate(grape_output)
    p_mag = plot_magnetization(iso_opt, t_c)
    p_ini = plot_magnetization(iso, t_c)
    p_field = plot_control_fields(fields_opt)
    p_init = plot_control_fields(field_init)
    display(p_mag)
    display(p_field)
end