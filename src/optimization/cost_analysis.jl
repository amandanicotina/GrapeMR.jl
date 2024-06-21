

labeloff  = ["T1 = 1000ms"]#, "White"];
targetoff = ["min"]
offset = 100           # [Hz]
offset_vals  = collect(-offset:2:offset)
n = 25
B1_vals = collect(range(0.7, stop=1.3, length=n))#length(offset_vals)))
spins_offset = GrapeMR.Spin(M0, [1.0], [0.6], offset_vals, [1.0], targetoff, labeloff)
spins_offset_B1 = GrapeMR.Spin(M0, [1.0], [0.6], offset_vals, B1_vals, targetoff, labeloff)
# spins_offset = GrapeMR.SteadyState(M0, [1.0], [0.6], offset_vals, [1.0], targetoff, labeloff, α, Δϕ, TR, TE)
# spins_offset_B1 = GrapeMR.SteadyState(M0, [1.0], [0.6], offset_vals, B1_vals, targetoff, labeloff, α, Δϕ, TR, TE)


# Calculate magnetization dynamics and cost function for each spin
isos = Isochromat[]
costs = []

for spin ∈ spins_offset
    mag = forward_propagation(grape_output.control_field, spin)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spin)
    push!(isos, iso)
    push!(costs, cost_function(iso, :saturation_contrast_My))
end
isos
p1 = plot_cost_offset(isos, :saturation_contrast_My)

# Extract B0 and B1 inhomogeneity values
B0_values = [spin.B0inho for spin in spins_offset_B1]
B1_values = [spin.B1inho for spin in spins_offset_B1]

# Prepare data for heatmap
cost_matrix = zeros(length(B1_values), length(B0_values))

@time for (i, b1) in enumerate(B1_values)
    for (j, b0) in enumerate(B0_values)
        s = GrapeMR.Spin(M0, [1.0], [0.6], b0, b1, targetoff, labeloff)
        mag = forward_propagation(grape_output.control_field, s[1])
        dyn = GrapeMR.Magnetization(mag)
        iso = Isochromat(dyn, s[1])
        cost_matrix[i, j] = cost_function(iso, :saturation_contrast_My)
    end
end

# Plot heatmap
using Plots

h1=heatmap(cost_matrix, xlabel="B0 Inhomogeneity", ylabel="B1 Inhomogeneity", title="Cost Function Map")
    xlims!(-100, 100)
    ylims!(0.7, 1.3)
display(p1) 
display(h1)