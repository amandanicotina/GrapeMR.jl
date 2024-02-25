using GrapeMR

##### INITIALIZATION #####
# SPINS #
M0 = [0.0; 0.0; 1.0];
T1 = [0.5];#, 0.7];
T2 = [0.2];#, 0.3];
target = ["max"];#, "min"];

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


##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, target, t_c, B1x, B1y, B0);


##### OPTIMIZE #####
opt_params = OptimizationParams(N, target_one_spin, [true true false]);
grape_output = grape(opt_params, field_init, spins; max_iter=10000, ϵ=1e-8); 
#grape_output_old = grape_optimize(opt_params, field_init, spins; max_iter=1000, ϵ=1e-3); 
##### PLOTS #####
PLOTS = false
if PLOTS
    for (idx, (iso_opt)) ∈ enumerate(grape_output.magnetization)
        p_mag = plot_magnetization(iso_opt, t_c)
        p_ini = plot_magnetization(iso_opt, t_c)
        p_field = plot_control_fields(grape_output.control_field)
        p_init = plot_control_fields(field_init)
        display(p_mag)
        #display(p_field)
    end
end

##### DEBUG #####
DEBUG = true
if DEBUG
    for (idx, spin) ∈ enumerate(spins)
        bx = plot(grape_output.cost_values[1,idx,:])
        by = plot(grape_output.cost_values[2,idx,:])
        bb = plot(grape_output.cost_values[1,idx,:])
             plot!(grape_output.cost_values[2,idx,:])
        #display(bx)
        #display(by)
        display(bb)
    ##### INITIAL MAGNETIZATION #####
        # mag = forward_propagation(field_init, spin);
        # dyn = Magnetization(mag)
        # iso = Isochromat(dyn, spin)
        # adj = backward_propagation(field_init, iso, grad_saturation_contrast)
        # push!(grape_output.isochromats, iso)
        # plot(adj')
        # ##### GRADIENTS #####
        # gx = eltype(Float64)[]
        # gy = eltype(Float64)[]
        # gx = gradient(adj, mag, Ix);
        # gy = gradient(adj, mag, Iy);
        # gxy = (gx, gy)
        # (bx, by) = update(field_init, gxy, 1e-6)
        # px = plot(gx')
        # py = plot(gy')
        # pbx = plot(bx')
        # pby = plot(by')
        # display(px)
        # display(py)
        # display(pbx)
        # display(pby)
    end
end
