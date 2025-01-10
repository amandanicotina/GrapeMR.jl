using GrapeMR

# Spin System
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
B0 = -200:5:200
target = ["max"]
label  = ["s3"]

T1 = [1.5]
T2 = [0.3]

α, Δϕ, TR, TE = 2π/9, 2π, 5e-3, 2.5e-3
spins = SteadyState(M0, T1, T2, B0, ΔB1, target, label, α, Δϕ, TR, TE) 

α = spins[1].α
offsets = collect(range(spins[1].B0inho, spins[end].B0inho, length(spins)))

plot_magnetization_trajectory(spins)


# Extract Mxy and Mz components
Mxy = steady_state_geometric.(spins)  # Transverse magnetization
Mz = steady_state_geometric_Mz.(spins)  # Longitudinal magnetization
angles = [π/6, π/4, π/3, π/2] 

plot_bssfp_magnetization(Mxy, Mz, offsets, α)

plot_ss_flip_angle_profile(spins, angles) 

plot_ss_offset_profile(spins) 



##### Dynamics #####
n = 100000
b = bSSFP_RF(n, 500, 2π/9, 5e-3)
plot_control_fields(b; unit="Hz")
iso = dynamics.(b, spins)

Mx = iso[1].magnetization.dynamics[2, :]
My = iso[1].magnetization.dynamics[3, :]
Mz = iso[1].magnetization.dynamics[4, :]
Mxy = sqrt.(Mx.^2 + My.^2)
ss_value = abs(steady_state(spins[1]))
t = range(0, b.t_control, length=n+1)

# Create plot
using Plots
p = plot(t, Mxy,
        xlabel="Time ms",
        ylabel="Transverse Magnetization |Mxy|",
        title="bSSFP Magnetization Evolution",
        label="Mxy",
        linewidth=2)
hline!([ss_value], 
linestyle=:dash, 
color=:red,
label="Steady State",
linewidth=1.5)
# Add annotation for steady state value
annotate!(0.3, ss_value-0.04, 
  text("SS = $(round(ss_value, digits=3))", 10, :red))

p = plot(t, Mz,
        xlabel="Time ms",
        ylabel="Longitudinal Magnetization |Mz|",
        title="bSSFP Longitudinal Magnetization Evolution",
        label="Mz",
        linewidth=2)





# Optimization Parameters
Tc = 1.0;
poly_start = 0.5;
poly_degree = 1.0;
max_iter = 2000;
opt_params   = OptimizationParams(poly_start, poly_degree, Int(ceil(max_iter)));

# Grape Parameters 
grape_params = GrapeParams(2000, GrapeMR.saturation_contrast_steady_state, Dict("B1x" => true, "B1y" => true, "Bz" => false));

# Parameters 
params = Parameters(grape_params, opt_params);

# Initial RF Pulse
B1ref = 1.0
control_field = spline_RF(grape_params.N, Tc, B1ref);

# Run Optimization
grape_output = @time grape(params, control_field, spins); 

# Plots
plot_cost_values(grape_output.cost_values, grape_params)
plot_magnetization_target(grape_output.isochromats)
plot_ss_offset_profile(grape_output.isochromats)
plot_magnetization_2D(grape_output.isochromats) 

plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)

plot_magnetization_target_3D(grape_output.isochromats[1])
plot_magnetization_target_3D(grape_output.isochromats[2])

plot_control_fields(grape_output.control_field; unit="Hz")


# Save Optimization Data
folder_path = "/Users/amandanicotina/Documents/PhD/Thesis/ChapterCodes/C4/"
save_grape_data(grape_output; folder_path)

