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
Tc, poly_start, poly_degree, max_iter = 0.836, 0.1, 1, 2000 #  bohb.minimizer # 0
opt_params   = OptimizationParams(poly_start, poly_degree,  Int(ceil(max_iter)))

# Grape Parameters 
grape_params = GrapeParams(1500, :saturation_contrast_Mx, [true true false])

# RF
B1ref = 1.0
B1x = spline_RF(grape_params.N, Tc)'
B1y = spline_RF(grape_params.N, Tc)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Optimization
grape_output = @time grape(opt_params, grape_params, control_field, spins);

# Plots
plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats);

## The Spin struct
# ''' julia
# struct Spin <: Spins
#   M_init::Vector{Float64}
#   T1::Float64
#   T2::Float64
#   #δ::Vector{Float64}
#   B0inho::Float64
#   B1inho::Float64
#   target::String
#   label::String
#   Nspins::Float64
# end
# '''