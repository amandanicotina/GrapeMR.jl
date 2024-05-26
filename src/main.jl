using GrapeMR
using BlochSim
using Hyperopt
using ParameterSchedulers

##### INITIALIZATION #####
# Parameters
N  = 1500     # Number of points
α  = 2π/9    # Flip angle in radians
Δϕ = π       # Phase cycling
TR = 4.28e-3 # Repetition time in seconds
TE = 2.19e-3 # Echo time in seconds
# Tc = 0.5     # Control time in seconds

bohb = @hyperopt for i = 50,
                     Tc = LinRange(0.1, 1, 10), 
                     poly_start = [1e-1, 1e-2, 1e-3, 1e-4], 
                     poly_degree = [1, 2, 3, 4];
                     sampler=Hyperband(R=50, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
    # BOHB
    # See: https://github.com/baggepinnen/Hyperopt.jl#bohb
    # And: https://arxiv.org/pdf/1807.01774
    # state is set by the BOHB algorithm and a KDE will estimate hyperparameters
    # that balance exploration and exploitation based on previous observations
    if state !== nothing
        Tc, poly_start, poly_degree = state
    end
    print(i, "\t", Tc, "\t", poly_start, "\t", poly_degree, "   \t")
    # RFs
    B1x = spline_RF(N, Tc)'
    B1y = spline_RF(N, Tc)'
    ΔB1, Bz = [1.0, 0.9, 0.8], zeros(1,N)

    # Spin parameters
    M0 = [0.0, 0.0, 1.0] # Initial magnetization vector
    T1 = [1.0, 0.08]     # Longitudinal relaxation time
    T2 = [0.6, 0.03]     # Transverse relaxation time
    B0 = 1/TR           # [Hz]
    B0_vals = [-B0, -B0/2, 0.0, B0/2, B0] 
    target  = ["min", "max"]
    label   = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

    # T1 = [1.830, 0.622, 2.430];
    # T2 = [0.184, 0.092, 0.132];
    # label = ["GD", "Yolk", "White"];
    # target = ["max", "min", "min"];

    # Spin and RF objects
    control_field = ControlField(B1x, B1y, 1.0, Bz, Tc)
    spins  = GrapeMR.Spin(M0, T1, T2, B0_vals, ΔB1, target, label)

    ##### OPTIMIZE #####
    max_iter     = 500
    lr_scheduler = Poly(start=poly_start, degree=poly_degree, max_iter=max_iter+1) 
    opt_params   = OptimizationParams(N, max_iter, :saturation_contrast, [true true false]);
    grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2);
    cost = grape_output.cost_values[end]
    @show cost
    cost, (Tc, poly_start, poly_degree)
end


# # RFs
# B1x = spline_RF(N, Tc)'
# B1y = spline_RF(N, Tc)'
# ΔB1, Bz = [1.0, 0.9, 0.8], zeros(1,N)

# # Spin parameters
# M0 = [0.0, 0.0, 1.0] # Initial magnetization vector
# T1 = [1.0, 0.08]     # Longitudinal relaxation time
# T2 = [0.6, 0.03]     # Transverse relaxation time
# B0 = 1/TR           # [Hz]
# B0_vals = [-B0, -B0/2, 0.0, B0/2, B0] 
# target  = ["min", "max"]
# label   = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

# # T1 = [1.830, 0.622, 2.430];
# # T2 = [0.184, 0.092, 0.132];
# # label = ["GD", "Yolk", "White"];
# # target = ["max", "min", "min"];

# # Spin and RF objects
# control_field = ControlField(B1x, B1y, 1.0, Bz, Tc)
# spins  = GrapeMR.Spin(M0, T1, T2, B0_vals, ΔB1, target, label)

# ##### OPTIMIZE #####
# max_iter     = 2000
# lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1) 
# opt_params   = OptimizationParams(N, max_iter, :saturation_contrast, [true true false]);
# grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 

##### PLOTS #####
#plot_magnetization_target_3d(grape_output.isochromats[2])
# plot_magnetization_2D(grape_output.isochromats) 
# plot_control_fields(grape_output.control_field) 
# plot_cost_values(grape_output.cost_values, opt_params)

##### SAVE DATA #####
# folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
# GrapeMR.save_grape_data(grape_output; folder_path)


# TODO Fix lr_scheduler
# TODO Check if lr_scheduler is actually working and ϵ is not just constant
# TODO Finish changing the cost_ss since now it's calculated inside SteadyState object
# TODO Think about the B0 issue since the SS changes for each different B0

