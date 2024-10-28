# # Tutorial
# The goal is to show how to implement the GrapeMR.jl package

# 1. *Define your physical system:* Spins, relaxation values, inhomogeneities, etc.
# 2. *Optimization parameters:* Scheduler parameters, max iterations, etc.
# 3. *Grape Parameters:* Time steps, cost function and which fields to optimize.
# 4. *Generate initial control field:* Spline 

using GrapeMR

# Physical system: Initial magnetization state, relaxation values in seconds, label are the spins names, 
# target: cost function dependent, B0 inhomonegenty in Herz. At the end, create the spin object with all spins
M0 = [0.0, 0.0, 1.0] 
T1 = [1.3, 0.73]
T2 = [2.0, 0.04]
label  = ["T1=$(T1[1]*1e3)ms", "T1=$(T1[2]*1e3)ms"]  
target = ["min", "max"]
B0 = 30.0
offset = collect(-B0/2:5:B0/2) 
ΔB1 = [1.0]
spins = Spin(M0, T1, T2, offset, ΔB1, target, label)

# Optimization Parameters:
max_iter = get(ENV, "DEV", "false") == "true" ? 1 : 2000  # we set max_iter to 1 if we're in development mode to build the docs faster
Tc, poly_start, poly_degree = 0.836, 0.1, 1
opt_params = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, :saturation_contrast_Mx, [true true false])

# RF
N = 1500
B1ref = 1.0
control_field = spline_RF(N, Tc, B1ref)

# Run Optimization
params = Parameters(grape_params, opt_params)
grape_output = grape(params, control_field, spins);

#md # # Plots
#md # ```@repl tutorial
#md # using GrapeMR, Plots; unicodeplots(); # change the backend so that plots go to stdout and can be rendered in CI/headless mode.
#md # default(show = false); #hide
#md # control_fields = plot_control_fields(grape_output.control_field);
#md # display(control_fields)
#md # ```

# The Spin struct:
[(k, v) for (k, v) in zip(fieldnames(GrapeMR.Spin), fieldtypes(GrapeMR.Spin))]
