# # Tutorial
# This tutorial demonstrates how to implement the GrapeMR.jl package by scripting a setup and directly calling its functions. 

# 1. *Define your physical system:* Spins, relaxation values, inhomogeneities, etc.
# 2. *Optimization parameters:* Scheduler parameters and maximum iterations.
# 3. *Grape Parameters:* Time steps, cost function and which fields to optimize.
# 4. *Generate initial control field:* Cubic spline is used by default.

using GrapeMR

# ## Saturation contrast without inhomogeneities

# ### Physical system

# Define the initial magnetization state, relaxation times (in seconds), and spin labels. 
# Specify `target` according to the desired outcome (cost function-dependent) and `offset` for $B_0$ inhomogeneities in Hertz. 
# Finally, create a [`Spin`](https://amandanicotina.github.io/GrapeMR.jl/dev/pages/api/#GrapeMR.Spin) object with all spins. Each unique inhomogeneity combination is treated as a separate spin.

M0 = [0.0, 0.0, 1.0]
T1 = [0.6, 0.1]
T2 = [0.3, 0.05]
ΔB1 = [1.0]
B0 = 0.0
offset = collect(-B0:1:B0)
label  = ["s1", "s2"]
target = ["min", "max"]
spins = Spin(M0, T1, T2, offset, ΔB1, target, label)

# ### Parameters

# Optimization Parameters
max_iter = get(ENV, "DEV", "false") == "true" ? 1 : 2000  # we set max_iter to 1 if we're in development mode to build the docs faster
poly_start = 0.5
poly_degree = 2
opt_params = OptimizationParams(poly_start, poly_degree, max_iter)

# Grape Parameters 
N = 2000
cost = GrapeMR.saturation_contrast
fields2opt = Dict("B1x" => true, "B1y" => true, "Bz" => false)
grape_params = GrapeParams(N, cost, fields2opt)

# `Parameter` struct
params = Parameters(grape_params, opt_params)

# ### Initial control field 

# Generate an initial control field using the spline function. Set the control time `control_time` and reference amplitude `B1ref`.

B1ref = 1.0
control_time = 0.5
control_field = spline_RF(N, control_time, B1ref)

# ### Run Optimization

# Run the optimization with the contructed `Spins`, configured `Parameters` and `ControlField`.

grape_output = grape(params, control_field, spins);

# ### Plot
#md # ```@repl tutorial
#md # using GrapeMR, Plots; unicodeplots(); # change the backend so that plots go to stdout and can be rendered in CI/headless mode.
#md # default(show = false); #hide
#md # control_fields = plot_control_fields(grape_output.control_field);
#md # display(control_fields)
#md # ```

# ## Magnetization saturation with $B_0$ inhomogeneity

# ### Physical system

# Optimizing pulses robust against field inhomogeneities. 
# Here $B_0$ is the field inhomogeneity in Hertz and $ΔB1$ corresponds to the RF field inhomogeneity in percentage.

M0 = [0.0, 0.0, 1.0]
T1 = [0.6]
T2 = [0.3]
ΔB1 = [1.0]
B0 = 15.0
label  = ["s1"]
target = ["saturation"]
offset = -B0:1:B0
spins = Spin(M0, T1, T2, offset, ΔB1, target, label)

# ### Parameters
# Optimization parameters
max_iter = get(ENV, "DEV", "false") == "true" ? 1 : 2000  # we set max_iter to 1 if we're in development mode to build the docs faster
poly_start = 0.5
poly_degree = 2
opt_params = OptimizationParams(poly_start, poly_degree, max_iter)

# Grape parameters
N = 2000
cost = GrapeMR.euclidean_norm
fields2opt = Dict("B1x" => true, "B1y" => true, "Bz" => false)
grape_params = GrapeParams(N, cost, fields2opt)

# `Parameter` struct
params = Parameters(grape_params, opt_params)

# ### Initial control fields

B1ref = 5.0
control_time = 0.5
control_field = spline_RF(N, control_time, B1ref)

# ### Run optimization
grape_output = grape(params, control_field, spins);

#md # ### Plots
#md # ```@repl tutorial
#md # using GrapeMR, Plots; unicodeplots(); # change the backend so that plots go to stdout and can be rendered in CI/headless mode.
#md # default(show = false); #hide
#md # control_fields = plot_control_fields(grape_output.control_field);
#md # display(control_fields)
#md # ```
