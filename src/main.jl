using GrapeMR
using BlochSim
using Hyperopt
using ParameterSchedulers

##### INITIALIZATION #####
bohb = @hyperopt for i = 50,
                     Tc = LinRange(0.05, 1.0, 20), 
                     poly_start = [1e-1, 1e-2], 
                     poly_degree = [1, 2, 3],
                     max_iter = range(2000, stop=8000, step=500);
                     sampler = Hyperband(R=50, η=3, inner=BOHB(dims=[Hyperopt.Continuous(), Hyperopt.Continuous(), Hyperopt.Continuous()])), 
    # BOHB
    # See: https://github.com/baggepinnen/Hyperopt.jl#bohb
    # And: https://arxiv.org/pdf/1807.01774
    # state is set by the BOHB algorithm and a KDE will estimate hyperparameters
    # that balance exploration and exploitation based on previous observations
    if state !== nothing
        Tc, poly_start, poly_degree, max_iter = state
    end
    print("\n", i, "\t", Tc, "\t", poly_start, "\t", poly_degree, "\t", max_iter, "\t")
    
    # RFs
    N   = 1000
    B1ref = 1.0
    B1x = spline_RF(N, Tc)'
    B1y = spline_RF(N, Tc)'
    Bz  = zeros(1,N)
    ΔB1 = [1.0]

    # Spin parameters
    M0 = [0.0, 0.0, 1.0] 
    T1 = [1.0, 0.08]#, 2.430];
    T2 = [0.6, 0.03]#, 0.132];
    label  = ["T1 = $(T1[1])", "T1 = $(T1[2])"]#, "White"];
    target = ["min", "max"]#, "min"];
    B0 = 25          
    offset = collect(-B0:1:B0) 
    
    # Spin and RF objects
    control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)
    spins  = GrapeMR.Spin(M0, T1, T2, offset, ΔB1, target, label)

    ##### OPTIMIZE #####
    lr_scheduler = Poly(start=poly_start, degree=poly_degree, max_iter=max_iter+1) 
    opt_params   = OptimizationParams(N, max_iter, :saturation_contrast, [true true false]);
    grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2);
    cost = grape_output.cost_values[end]
    @show cost
    cost, (Tc, poly_start, poly_degree)
end

# use bohb configuration parameters
Tc, poly_start, poly_degree, max_iter = bohb.minimizer
# Tc, poly_start, poly_degree, max_iter = 0.8, 1e-2, 3, 5000

# RFs
N   = 1000
B1ref = 1.0
B1x = spline_RF(N, Tc)'
B1y = spline_RF(N, Tc)'
Bz  = zeros(1,N)
ΔB1 = [1.0]

# Spin parameters
M0 = [0.0, 0.0, 1.0] 
T1 = [1.0, 0.08]#, 2.430];
T2 = [0.6, 0.03]#, 0.132];
label  = ["T1 = $(T1[1])", "T1 = $(T1[2])"]#, "White"];
target = ["min", "max"]#, "min"];
B0 = 25          
offset = collect(-B0:1:B0) 

# Spin and RF objects
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)
spins  = GrapeMR.Spin(M0, T1, T2, offset, ΔB1, target, label)

# Spin and RF objects
control_field = ControlField(B1x, B1y, B1ref, Bz, Tc)
spins  = GrapeMR.Spin(M0, T1, T2, offset, ΔB1, target, label)

##### OPTIMIZE #####
lr_scheduler = Poly(start=poly_start, degree=poly_degree, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, :saturation_contrast, [true true false])
grape_output = @time grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2);

#### PLOTS #####
#plot_magnetization_target_3D(grape_output.isochromats[2])
plot_magnetization_2D(grape_output.isochromats) 
plot_control_fields(grape_output.control_field) 
plot_cost_offset(grape_output.isochromats, :saturation_contrast)
plot_cost_values(grape_output.cost_values, opt_params)
plot_magnetization_time(grape_output.isochromats[1], Tc)

(pCost, pOrder) = bohb_params(bohb)
display(pCost)
display(pOrder)

# Save data
# folder_path = "/Users/amandanicotina/Documents/Documents/ProgressReports/ResultsGrapeMR"
# GrapeMR.save_grape_data(grape_output; folder_path)
# using Plots
# pMag = plot(titlefontsize = 14)

labeloff  = ["T1 = 80ms"]#, "White"];
targetoff = ["min"]
offset = 25           # [Hz]
offset_vals  = collect(-offset:0.1:offset) 
B1_vals = collect(range(0.7, stop=1.3, length=length(offset_vals)))
spins_offset = GrapeMR.Spin(M0, [0.08], [0.03], offset_vals, B1_vals, targetoff, labeloff)

# Calculate magnetization dynamics and cost function for each spin
isos = []
costs = []

for spin ∈ spins_offset
    mag = forward_propagation(grape_output.control_field, spin)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spin)
    push!(isos, iso)
    push!(costs, cost_function(iso, :saturation_contrast))
end

# Extract B0 and B1 inhomogeneity values
B0_values = [spin.B0inho for spin in spins_offset]
B1_values = [spin.B1inho for spin in spins_offset]

# Prepare data for heatmap
# Assuming B0 and B1 values are regularly spaced, otherwise you need interpolation
B0_unique = unique(B0_values)
B1_unique = unique(B1_values)
cost_matrix = zeros(length(B1_unique), length(B0_unique))

for (i, b1) in enumerate(B1_unique)
    for (j, b0) in enumerate(B0_unique)
        s = GrapeMR.Spin(M0, [1.0], [0.6], b0, b1, targetoff, labeloff)
        mag = forward_propagation(grape_output.control_field, s[1])
        dyn = GrapeMR.Magnetization(mag)
        iso = Isochromat(dyn, s[1])
        cost_matrix[i, j] = cost_function(iso, :saturation_contrast)
    end
end

# Plot heatmap
heatmap(B0_unique, B1_unique, cost_matrix, xlabel="B0 Inhomogeneity", ylabel="B1 Inhomogeneity", title="Cost Function Map")

##### Phase testing 
# ux = grape_output.control_field.B1x
# uy = grape_output.control_field.B1y
# ψ  = π
# Ux = vec(ux*cos(ψ) .+ uy*sin(ψ)) 
# Uy = vec(uy*cos(ψ) .- ux*sin(ψ))

# tc = grape_output.control_field.t_control
# t = range(0.0, tc, length = length(ux))
# p_Bx = plot(t, Ux, linewidth=2, label=false, ylabel="B1x", title="Control Fields", titlefontsize=15)
# p_By = plot(t, Uy, linewidth=2, label=false, ylabel="B1y", xlabel="t [s]")
# # xticks_values = [-π, -π/2, 0, π/2, π]
# # xticks_labels = ["-π", "-π/2", "0", "π/2", "π"]
# # p_By = plot!(p_By, xticks=(xticks_values, xticks_labels))

# p = plot(p_Bx, p_By, layout = (2,1))

# plot(vec(Ux), vec(Uy))