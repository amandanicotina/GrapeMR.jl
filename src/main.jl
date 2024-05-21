using GrapeMR
using BlochSim
using ParameterSchedulers

##### INITIALIZATION #####
# Parameters
N  = 1000     # Number of points
α  = 2π/9    # Flip angle in radians
Δϕ = π       # Phase cycling
TR = 4.28e-3 # Repetition time in seconds
TE = 2.19e-3 # Echo time in seconds       
Tc = 0.5     # Control time in seconds

# RFs
B1x = spline_RF(N, Tc)'
B1y = spline_RF(N, Tc)'
ΔB1, Bz = [1.0], zeros(1,N)

# Spin parameters
M0 = [0.0, 0.0, 1.0] # Initial magnetization vector
T1 = [1.0, 0.08]     # Longitudinal relaxation time
T2 = [0.6, 0.03]     # Transverse relaxation time
B0 = [100.0]           # [Hz]
B0_vals = range(-B0, B0, length=N) 
target  = ["max", "min"]
label   = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

# T1 = [1.830, 0.622, 2.430];
# T2 = [0.184, 0.092, 0.132];
# label = ["GD", "Yolk", "White"];
# target = ["max", "min", "min"];

# Spin and RF objects
control_field = ControlField(B1x, B1y, 1.0, Bz, Tc)
spins  = GrapeMR.SteadyState(M0, T1, T2, B0, ΔB1, target, label, α, Δϕ, TR, TE)

##### OPTIMIZE #####
max_iter     = 3000
lr_scheduler = Poly(start=1e-1, degree=1, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, :target_steady_state, [true true false]);
grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 

##### PLOTS #####
plot_magnetization_target_3d(grape_output.isochromats[2])
plot_magnetization_target(grape_output.isochromats)
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

##### SAVE DATA #####
folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
GrapeMR.save_grape_data(grape_output; folder_path)


# TODO Fix lr_scheduler
# TODO Check if lr_scheduler is actually working and ϵ is not just constant
# TODO Finish changing the cost_ss since now it's calculated inside SteadyState object
# TODO Think about the B0 issue since the SS changes for each different B0

