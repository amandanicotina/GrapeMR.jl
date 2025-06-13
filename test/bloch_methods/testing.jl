using GrapeMR
using Plots

# α = π / 2
# t_c = 1e-4
# B1 = α / (2π * t_c)  # Negated to rotate into +My

B1 = 1.0
t_c = 0.5

cf_90x = generate_control_field(:hard; t_c=t_c, B1ref=B1, normalize=true)

spin_rot = Spin([0.0, 0.0, 1.0], 1e8, 1e8, 0.0, 1.0, "mz", "test", 1)
M_rot = zeros(4, size(cf_90x.B1x, 2) + 1)
forward_propagation!(M_rot, cf_90x, spin_rot)

time = LinRange(0.0, t_c, length(cf_90x.B1x))
plot(layout = (2,1))
plot!(time, cf_90x.B1x', ylabel = "B1x [Hz]", label = false, subplot = 1)
plot!(time, cf_90x.B1y', ylabel = "B1y [Hz]", xlabel = "time [sec]", label = false, subplot = 2)

time1 = LinRange(0.0, t_c, length(cf_90x.B1x) + 1)
plot(time1, M_rot[2:end,:]')


# Mx = M_rot[2,end]
# My = M_rot[3,end]
# Mz = M_rot[4,end]
# sqrt(Mx^2 + My^2 + Mz^2 + 1e-12)

# Spin Parameters
M0 = [0.0, 0.0, 1.0]
ΔB1 = [1.0]
offsets = -15:15:0
T1 = [1e8] #[1/31.3436]
T2 = [1e8] #[1/37.6471]
label = ["C1"]
target = ["min"]
spins_90x = generate_spins(M0, T1, T2, offsets, ΔB1, target, label)

B1 = 1.0
t_c = 0.5
cf_90x = generate_control_field(:hard; t_c=t_c, B1ref=B1, normalize=true)
magnetization = []
for spin in spins_90x
    M_rot = zeros(4, size(cf_90x.B1x, 2) + 1)
    forward_propagation!(M_rot, cf_90x, spin)
    push!(magnetization, M_rot)
end

time1 = LinRange(0.0, t_c, length(cf_90x.B1x) + 1)
p = plot(xlabel = "time", ylabel = "Magnetization")
for n in 1:length(spins_90x)
    plot!(p, time1, magnetization[n][4,:])
    scatter!(p, [time1[end]], [magnetization[n][4,end]], label =false)
end
display(p)

n = length(cf_90x.B1x)
magnetization = zeros(Float64, 4, n + 1)
adjoint = zeros(Float64, 4, n + 1)
forward_propagation!(magnetization, cf_90x, spins_90x[1])
dyn = Magnetization(magnetization)
iso = Isochromat(dyn, spins_90x[1])

cost, adj_init = grape_params.cost_function(iso)
backward_propagation!(adjoint, cf_90x, iso, adj_init)
∇x = zeros(Float64, 1, n)
∇x .+= gradient(adjoint, magnetization, Ix)
∇y = zeros(Float64, 1, n)
∇y .+= gradient(adjoint, magnetization, Iy)
ϵ=1e-3

update!(cf_90x, (∇x, ∇y), ϵ)

@code_warntype bloch_matrix(cf_90x.B1x[1], cf_90x.B1y[1], cf_90x.Bz[1], T1[1], T2[1])