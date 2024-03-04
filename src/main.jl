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
# N   = 600;
# t_c = 1.0;
# B1  = 10.0;
# B1x, B1y = B1*ones(Float64, 1, N), 5.0*ones(Float64, 1, N)
# B0, ΔB0  = zeros(1, N), zeros(1, N);

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
B0    = zeros(1,N);#[-150.0 0.0 150.0]; #Hz
B1x   = ((flip_x.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
B1y   = ((flip_y.*sinc.(y))./2π)';
#B1y   = zeros(1, N);
#plot(B1x')

##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, target, label, t_c, B1x, B1y, B0);


##### OPTIMIZE #####
opt_params = OptimizationParams(N, cost_saturation_contrast, [true true false]);
grape_output = @time grape(opt_params, field_init, spins; max_iter=25000, ϵ=1e-3); 


##### PLOTS #####
PLOTS = false
if PLOTS
    for (idx, (iso_opt)) ∈ enumerate(grape_output.magnetization.dynamics)
        p_mag = plot_magnetization(iso_opt, t_c)
        p_ini = plot_magnetization(iso_opt, t_c)
        p_field = plot_control_fields(grape_output.control_field)
        p_init = plot_control_fields(field_init)
        display(p_mag)
        #display(p_field)
    end
end

##### DEBUG #####
DEBUG = false
if DEBUG
    for (idx, spin) ∈ enumerate(spins)
        ##### GRADIENTS #####
        gx = eltype(Float64)[]
        gy = eltype(Float64)[]
        gx = gradient(adj, mag, Ix);
        gy = gradient(adj, mag, Iy);
        gxy = (gx, gy)
        (bx, by) = update(field_init, gxy, 1e-4)
        px = plot(gx')
        py = plot(gy')
        pbx = plot(bx')
        pby = plot(by')
        display(px)
        #display(py)
        #display(pbx)
        #display(pby)
    end
end

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