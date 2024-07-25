# using Distributed
# const num_processes::Int = 6;
# addprocs(num_processes)
# println(nprocs())
# println(workers())

# @everywhere begin
using GrapeMR
using Plots
using JLD2
using Hyperopt

M0 = [0.0, 0.0, 1.0]
T1 = [1.830] # collect(range(10, stop = 250, length=100).*1e-3)#, 1.0]
T2 = [0.184] # collect(range(10, stop = 250, length=100).*1e-3)#, 0.6]
label  = fill("-", length(T1)) 
target = fill("min", length(T1))
ΔB1 = [1.0]
α, Δϕ, TR, TE = 2π/9, π, 5.6e-3, 2.8e-3
B0 = round(1/TR)
offset = collect(-B0:5:B0) 
spins = GrapeMR.SteadyState(M0, T1, T2, offset, ΔB1, target, label, α, Δϕ, TR, TE)
ss = steady_state_matrix.(spins)
Mx = getproperty.(ss, :x)
My = getproperty.(ss, :y)
Mxy = sqrt.(Mx.^2 + My.^2)
# spins = GrapeMR.Spin(M0, T1, T2, [0.0], ΔB1, target, label)
plot(offset, Mxy)

# Grape Parameters 
grape_params = GrapeParams(1000, :target_steady_state, [true true false])
# end

# const wandb_project::String = "GrapeMR"
# logger::WandbLogger = Wandb.WandbLogger(; project = wandb_project, name = nothing)
# bohb = @time hyperoptimization(spins, grape_params, LinRange(0.05, 1.5, 100), 5000, i=100)
# # rand_hopt = random_sample(spins, grape_params, LinRange(0.05, 1.0, 20), range(1000, stop=5000, step=100), i=20)
# # close(logger)

# p1, p2 = plot_bohb(bohb)
# display(p1)

# Optimization Parameters
# Tc, poly_start, poly_degree, max_iter = 0.5, 1e-1, 2, 2000 # bohb.minimizer # 
# cost = 0.006355491696394037
# cost = 0.006639151222756831
Tc = 1.5
poly_start = 0.1
poly_degree = 1.0
max_iter = 5000
opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)))

# RF
B1ref = 1.0
B1x = spline_RF(grape_params.N, Tc)'
B1y = spline_RF(grape_params.N, Tc)'
Bz  = zeros(1, grape_params.N)
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)

# Run Grape
grape_output_1 = @time grape(opt_params, grape_params, control_field, spins);
grape_output_2 = @time grape(opt_params, grape_params, control_field, spins);
grape_output_3 = @time grape(opt_params, grape_params, control_field, spins);

grape_output_1.cost_values[end]
plot_cost_values(grape_output_1.cost_values, grape_params)
plot_magnetization_2D(grape_output_1.isochromats)
plot_magnetization_control_field(grape_output_1.control_field, grape_output_euromar.isochromats)
plot_Mtrans_offset_ss(grape_output_1.isochromats)


JLD2.@save "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR/2024-07-25/bssfp_opt.jl" grape_output_3 opt_params grape_params


# Dynamics outside the optimized range
M0 = [0.0, 0.0, 1.0]
T1 = collect(range(1, stop = 500, length=500).*1e-3)#[1.830]#, 1.0]
T2 = collect(range(1, stop = 500, length=500).*1e-3)#[0.184]#, 0.6]
label  = fill("-", length(T1)) 
target = fill("min", length(T1))
B0 = 30.0
offset = collect(-B0/2:5:B0/2) 
ΔB1 = [1.0]

spins_outside = GrapeMR.Spin(M0, T1, T2, [0.0], ΔB1, target, label)
cf = grape_output_euromar.control_field
Base.broadcastable(cf::ControlField) = Ref(cf)
mag_test = forward_propagation.(cf, spins_outside)
dyn_test = GrapeMR.Magnetization.(mag_test)
iso_test = Isochromat.(dyn_test, spins_outside)
plot_magnetization_2D(iso_test)



# cost = 0.006355491696394037
# Tc = 1.2217171717171718
# poly_start = 0.1
# poly_degree = 1.0
# max_iter = 5000

cost = 0.006639151222756831
Tc = 1.5
poly_start = 0.1
poly_degree = 1.0
max_iter = 5000

# cost = 0.16551798814860533
# Tc = 0.40151515151515155
# poly_start = 0.1
# poly_degree = 2.0
# max_iter = 5000

# cost = 0.013969845439110222
# Tc = 0.8409090909090908
# poly_start = 0.1
# poly_degree = 2.0
# max_iter = 5000