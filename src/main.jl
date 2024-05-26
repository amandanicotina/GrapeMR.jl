using GrapeMR
using BlochSim
using ParameterSchedulers

##### INITIALIZATION #####
# Parameters
N  = 1000     # Number of points
Tc = 1.0     # Control time in seconds

# RFs
B1x = spline_RF(N, Tc)'
B1y = spline_RF(N, Tc)'
ΔB1, Bz = [1.0], zeros(1,N)

# Spin parameters
# T1 = [1.0, 0.08]     # Longitudinal relaxation time
# T2 = [0.6, 0.03]     # Transverse relaxation time
# target  = ["min", "max"]
# label   = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

M0 = [0.0, 0.0, 1.0] 
T1 = [0.8, 0.4]#, 2.430]
T2 = [0.6, 0.2]#, 0.132]
label = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]#, "White"]
target = ["min", "max"]#, "min"]
B0 = 50.0
B0_vals = [-B0, -B0/2, 0.0, B0/2, B0] 

# Spin and RF objects
control_field = ControlField(B1x, B1y, 1.0, Bz, Tc)
spins  = GrapeMR.Spin(M0, T1, T2, B0_vals, ΔB1, target, label)

##### OPTIMIZE #####
max_iter     = 15000
lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, :saturation_contrast, [true true false]);
grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 

##### PLOTS #####
#plot_magnetization_target_3d(grape_output.isochromats[2])
plot_magnetization_2D(grape_output.isochromats) 
plot_control_fields(grape_output.control_field) 
plot_cost_values(grape_output.cost_values, opt_params)

##### SAVE DATA #####
folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
GrapeMR.save_grape_data(grape_output; folder_path)


# TODO Fix lr_scheduler
# TODO Check if lr_scheduler is actually working and ϵ is not just constant
# TODO Finish changing the cost_ss since now it's calculated inside SteadyState object
# TODO Think about the B0 issue since the SS changes for each different B0

