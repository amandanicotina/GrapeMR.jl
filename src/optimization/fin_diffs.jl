using GrapeMR
using Plots

# Parameters
t_c = 0.1                    # Total duration (s)
B1ref = 4.0                  # Tesla reference field
δ = 1e-6                     # Finite difference step (in Tesla)
cf = generate_control_field(:spline; t_c=t_c, B1ref=B1ref)
N = size(cf.B1x, 2)

# Spin setup
m_init = [0.0, 0.0, 1.0]
spin = Spin(m_init, 1.0, 0.2, 0.0, 1.0, "origin", "test", 1)

# Forward propagation
M = zeros(4, N + 1)
forward_propagation!(M, cf, spin)
iso = Isochromat(GrapeMR.Magnetization(M), spin)
val, cost_grad = euclidean_norm(iso)

# Backward propagation
χ = zeros(4, N + 1)
backward_propagation!(χ, cf, iso, cost_grad)

# GRAPE (adjoint) gradients
grad_x = zeros(1, N)
grad_y = zeros(1, N)
gradient!(grad_x, χ, M, Ix)
gradient!(grad_y, χ, M, Iy)

# Finite difference gradients (normalized steps)
grad_fd_x = zeros(1, N)
grad_fd_y = zeros(1, N)

δ_norm = δ / cf.B1_ref  # Normalized delta
for i in 1:N
    # B1x perturbed +δ
    cf_p = deepcopy(cf)
    cf_p.B1x[1, i] += δ_norm
    M_p = zeros(4, N + 1)
    forward_propagation!(M_p, cf_p, spin)
    val_p, _ = euclidean_norm(Isochromat(Magnetization(M_p), spin))

    # B1x perturbed -δ
    cf_m = deepcopy(cf)
    cf_m.B1x[1, i] -= δ_norm
    M_m = zeros(4, N + 1)
    forward_propagation!(M_m, cf_m, spin)
    val_m, _ = euclidean_norm(Isochromat(Magnetization(M_m), spin))

    grad_fd_x[1, i] = (val_p - val_m) / (2δ)
end

for i in 1:N
    # B1y perturbed +δ
    cf_p = deepcopy(cf)
    cf_p.B1y[1, i] += δ_norm
    M_p = zeros(4, N + 1)
    forward_propagation!(M_p, cf_p, spin)
    val_p, _ = euclidean_norm(Isochromat(Magnetization(M_p), spin))

    # B1y perturbed -δ
    cf_m = deepcopy(cf)
    cf_m.B1y[1, i] -= δ_norm
    M_m = zeros(4, N + 1)
    forward_propagation!(M_m, cf_m, spin)
    val_m, _ = euclidean_norm(Isochromat(Magnetization(M_m), spin))

    grad_fd_y[1, i] = (val_p - val_m) / (2δ)
end

# Diagnostics
Δx = grad_x .- grad_fd_x
Δy = grad_y .- grad_fd_y
rel_err_x = maximum(abs.(Δx)) / maximum(abs.(grad_x))
rel_err_y = maximum(abs.(Δy)) / maximum(abs.(grad_y))

println("\n✅ Normalized Finite Difference Gradient Check")
println("Max ∇ₓ GRAPE     = ", maximum(abs.(grad_x)))
println("Max ∇ₓ FinDiff   = ", maximum(abs.(grad_fd_x)))
println("Relative Error x = ", rel_err_x)

println("Max ∇ᵧ GRAPE     = ", maximum(abs.(grad_y)))
println("Max ∇ᵧ FinDiff   = ", maximum(abs.(grad_fd_y)))
println("Relative Error y = ", rel_err_y)

# Plot results
time = range(0, 1, length=N)

plot(layout = (2, 2), size=(700, 500))

plot!(time, grad_x[1, :], lw=2, label="GRAPE ∇ₓ", title="∇ₓ Comparison", xlabel="Normalized Time", ylabel="Gradient", subplot=1)
plot!(time, grad_fd_x[1, :], lw=2, linestyle=:dash, label="FiniteDiff ∇ₓ", subplot=2)

plot!(time, grad_y[1, :], lw=2, label="GRAPE ∇ᵧ", title="∇ᵧ Comparison", xlabel="Normalized Time", ylabel="Gradient", subplot=3)
plot!(time, grad_fd_y[1, :], lw=2, linestyle=:dash, label="FiniteDiff ∇ᵧ", subplot=4)


using GrapeMR
# using Plots

# # Parameters
# t_c = 0.1
# B1ref = 4.0          # Reference B1 amplitude (Hz)
# δ = 1e-5             # Perturbation for finite difference (unitless)

# # Generate normalized control field
# cf = generate_control_field(:gaussian; t_c=t_c, B1ref=B1ref)
# @assert cf isa NormalizedControlField "Expected a NormalizedControlField"
# @assert isapprox(cf.B1_ref, 1.0; atol=1e-12) "Control field is not normalized"

# # Check field scale
# println("→ Max B1x: ", maximum(abs, cf.B1x))
# println("→ Max B1y: ", maximum(abs, cf.B1y))
# plot(cf.B1x', label="B1x (normalized)", xlabel="Time step", ylabel="Amplitude")

# # Spin system setup
# m_init = [0.0, 0.0, 1.0]
# spin = Spin(m_init, 1.0, 0.2, 0.0, 1.0, "origin", "test", 1)

# # Forward propagation
# N = length(cf.B1x)
# M = zeros(4, N + 1)
# forward_propagation!(M, cf, spin)
# iso = Isochromat(GrapeMR.Magnetization(M), spin)
# val, cost_grad = euclidean_norm(iso)

# # Backward propagation
# χ = zeros(4, N + 1)
# backward_propagation!(χ, cf, iso, cost_grad)

# # GRAPE gradients
# grad_x = zeros(1, N)
# grad_y = zeros(1, N)
# gradient!(grad_x, χ, M, Ix)
# gradient!(grad_y, χ, M, Iy)

# # Finite difference gradients
# grad_fd_x = zeros(1, N)
# grad_fd_y = zeros(1, N)

# for i in 1:N
#     # Perturb B1x directly in normalized units
#     B1x_p = copy(cf.B1x); B1x_p[i] += δ
#     B1x_m = copy(cf.B1x); B1x_m[i] -= δ

#     cf_p = NormalizedControlField(B1x_p, cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)
#     cf_m = NormalizedControlField(B1x_m, cf.B1y, cf.B1_ref, cf.Bz, cf.t_control)

#     M_p = zeros(4, N + 1); forward_propagation!(M_p, cf_p, spin)
#     M_m = zeros(4, N + 1); forward_propagation!(M_m, cf_m, spin)

#     iso_p = Isochromat(GrapeMR.Magnetization(M_p), spin)
#     iso_m = Isochromat(GrapeMR.Magnetization(M_m), spin)

#     val_p, _ = euclidean_norm(iso_p)
#     val_m, _ = euclidean_norm(iso_m)

#     grad_fd_x[1, i] = (val_p - val_m) / (2δ)
# end



# for i in 1:N
#     cf_p = deepcopy(cf); cf_p.B1y[i] += δ 
#     M_p = zeros(4, N + 1)
#     forward_propagation!(M_p, cf_p, spin)
#     iso_p = Isochromat(GrapeMR.Magnetization(M_p), spin)
#     val_p, _ = euclidean_norm(iso_p)

#     cf_m = deepcopy(cf); cf_m.B1y[i] -= δ
#     M_m = zeros(4, N + 1)
#     forward_propagation!(M_m, cf_m, spin)
#     iso_m = Isochromat(GrapeMR.Magnetization(M_m), spin)
#     val_m, _ = euclidean_norm(iso_m)

#     grad_fd_y[1, i] = (val_p - val_m) / (2δ)
# end

# # Diagnostics
# max_Δx = maximum(abs.(grad_x .- grad_fd_x))
# rel_err_x = max_Δx / maximum(abs.(grad_x))

# max_Δy = maximum(abs.(grad_y .- grad_fd_y))
# rel_err_y = max_Δy / maximum(abs.(grad_y))

# println("\n✅ Normalization Check Passed")
# println("Max ∇x GRAPE     = ", maximum(abs.(grad_x)))
# println("Max ∇x FinDiff   = ", maximum(abs.(grad_fd_x)))
# println("Relative Error x = ", rel_err_x)
# println()
# println("Max ∇y GRAPE     = ", maximum(abs.(grad_y)))
# println("Max ∇y FinDiff   = ", maximum(abs.(grad_fd_y)))
# println("Relative Error y = ", rel_err_y)

# # Optional: Plot comparison
# plot(grad_x[1, :], label="∇x GRAPE")
# plot(grad_fd_x[1, :], label="∇x FiniteDiff", ls=:dash)

# plot(grad_y[1, :], label="∇y GRAPE")
# plot!(grad_fd_y[1, :], label="∇y FiniteDiff", ls=:dash)

# plot(grad_y[:], label="Adjoint ∇B₁y")
# plot!(grad_fd_y[:], label="FD ∇B₁y", linestyle=:dash)
# title!("Gradient Check: B1y")


# @testset "Analytical Gradient vs Finite Differences" begin
#     # Setup: control field
#     t_c = 0.01
#     B1ref = 1.0
#     cf = generate_control_field(:spline; t_c=t_c, B1ref=B1ref)

#     # Spin setup
#     m_init = [0.0, 0.0, 1.0]
#     spin = Spin(m_init, 1.0, 0.1, 0.0, 1.0, "mz", "spin", 1)
#     iso = Isochromat(GrapeMR.forward_propagation(cf, spin), spin)

#     # Gradient direction (e.g., maximize Mz)
#     grad = [0.0, 0.0, 0.0, 1.0]

#     # Compute analytical gradient
#     χ = zeros(4, size(cf.B1x, 2) + 1)
#     backward_propagation!(χ, cf, iso, grad)
#     ∇an = analytical_gradient(cf, iso, χ)  # Tuple (∇x, ∇y)

#     # Finite differences
#     ∇fd = finite_difference_field(cf, spin, grad; ϵ=1e-6)  # Tuple (∇x_fd, ∇y_fd)

#     # Compare gradients
#     @test isapprox(∇an[1], ∇fd[1]; rtol=1e-2)  # x component
#     @test isapprox(∇an[2], ∇fd[2]; rtol=1e-2)  # y component
# end
