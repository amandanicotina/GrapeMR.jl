using GrapeMR
using BlochSim
using Plots

# Parameters

# RFs
N  = 401                 # Number of points
α  = π / 8               # Flip angle in radians
Δϕ = π                  # Phase cycling
TR = 5e-3               # Repetition time in seconds

# Spin parameters
M0 = [0.0, 0.0, 1.0]    # Initial magnetization vector
T1 = [1.0]              # Longitudinal relaxation time
T2 = [0.8]              # Transverse relaxation time
B0 = range(-2/TR, stop=2/TR, length=N) |> collect # B0 field range

# Target and label for simulation
target = ["max"]
label  = ["T1-$(Int(T1[1]*1e3))ms"]

# Create spins using SteadyState method from GrapeMR
spins = GrapeMR.SteadyState(M0, T1, T2, B0, 1.0, target, label, α, Δϕ, TR, TR/2)

# Preallocate arrays for results
sig_abs   = Vector{Float64}(undef, N)
sig_phase = Vector{Float64}(undef, N)

sig_geo_abs = Vector{Float64}(undef, N)
sig_geo_Mz  = Vector{Float64}(undef, N)

sig_mat_abs = Vector{Float64}(undef, N)
sig_mat_phase = Vector{Float64}(undef, N)
mag_mat = Vector{BlochSim.Magnetization}(undef, N)

# Simulation
for (i, spin) ∈ enumerate(spins)
    ss = steady_state(spin)
    sig_abs[i]   = abs(ss)
    sig_phase[i] = angle(ss)
    
    ss_geo = steady_state_geometric(spin)
    sig_geo_abs[i] = abs(ss_geo)

    # ss_geo_Mz = steady_state_geometric_Mz(spin)
    # sig_geo_Mz[i] = abs(ss_geo_Mz)
    
    ss_mat, M_mat = steady_state_matrix(spin)
    sig_mat_abs[i] = abs(ss_mat)
    sig_mat_phase[i] = angle(ss_mat)
    mag_mat[i] = M_mat
end

# Plotting results
pPhase = plot(B0, sig_phase, label = "Manual", ylabel = "Phase [rad]")
         plot!(pPhase, B0, sig_mat_phase, label = "Matrix")

pMag = plot(B0, sig_abs, label = "Manual", xlabel = "Offset [Hz]", ylabel = "Magnitude")
       plot!(pMag, B0, sig_geo_abs, label = "Geometric")
       plot!(pMag, B0, sig_mat_abs, label = "Matrix")

pSig = plot(pMag, pPhase; layout = (2,1))


plot(B0, getproperty.(mag_mat, :z), label = "Matrix")
plot!(B0, Mzgeo, label = "Geometric", xlabel = "Offset [Hz]", ylabel="Mz")

display(pSig)

# Tests 