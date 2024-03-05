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
t_c = 0.6;

time = range(0.0, t_c, N);
t    = time .- t_c/2;
rotx = rad2deg(αx)/360;
roty = rad2deg(αy)/360;
flip_x = rotx/diff(t)[1];
flip_y = roty/diff(t)[1];
BW_Hz = 300.0;
x     = BW_Hz.*t;
y     = BW_Hz.*t;
#B0    = [0.0];
B0    = [-150.0, -100.0, -50, 0.0, 50.0, 100.0, 150.0]; # [Hz]
Bz    = zeros(1,N);
B1x   = ((flip_x.*sinc.(x))./2π)'; # sinc(x) = sin(πx)/(πx)
B1y   = ((flip_y.*sinc.(y))./2π)';

##### NORMALIZE #####
(spins, field_init) = normalization(M0, T1, T2, B0, target, label, t_c, B1x, B1y, Bz);
plot_control_fields(field_init)

##### OPTIMIZE #####
opt_params = OptimizationParams(N, cost_saturation_contrast, [true true false]);
grape_output = @time grape(opt_params, field_init, spins; max_iter=15000, ϵ=1e-2); 

##### PLOTS #####
plot_magnetization_Mz_Mxy(grape_output.isochromats)
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

##### Weigths & Biases #####
if wandb == false
    using Wandb, Dates, Logging

    # Start a new run, tracking hyperparameters in config
    lg = WandbLogger(project = "GrapeMR.jl", name = "grapemr-demo-$(now())",
                    config = Dict("learning_rate" => 0.01, 
                                "dropout" => 0.2,
                                "architecture" => "CNN", 
                                "dataset" => "CIFAR-100"))

    # Use LoggingExtras.jl to log to multiple loggers together
    global_logger(lg)

    # Simulating the training or evaluation loop
    for x in 1:50
        acc = log(1 + x + rand() * get_config(lg, "learning_rate") + rand() + 
                    get_config(lg, "dropout"))
                    
        loss = 10 - log(1 + x + rand() + x * get_config(lg, "learning_rate") + rand() +
                        get_config(lg, "dropout"))
        # Log metrics from your script to W&B
        @info "metrics" accuracy=acc loss=loss
    end

    # Finish the run
    close(lg)
end