using GrapeMR
using BlochSim
using Plots

# Parameters
# RFs
N  = 501     # Number of points
α  = 2π/9    # Flip angle in radians
Δϕ = π       # Phase cycling
TR = 4.28e-3 # Repetition time in seconds
TE = 2.19e-3 # Echo time in seconds

# Spin parameters
M0 = [0.0, 0.0, 1.0]
T1 = [1.83]         
T2 = [0.18]          
B0 = 100       
B0_vals = range(-B0, B0, length=N) 

# Target and label for simulation
target = ["max"]
label  = ["T1-$(Int(T1[1]*1e3))ms"]

# Create spins using SteadyState method from GrapeMR
spins = GrapeMR.SteadyState(M0, T1, T2, B0_vals, 1.0, target, label, α, Δϕ, TR, TE)

# Preallocate arrays for results
ss_vec_Mz    = Vector{Float64}(undef, N)
ss_vec_abs   = Vector{Float64}(undef, N)
ss_vec_phase = Vector{Float64}(undef, N)

sig_abs   = Vector{Float64}(undef, N)
sig_phase = Vector{Float64}(undef, N)

sig_geo_abs = Vector{Float64}(undef, N)
sig_geo_Mz  = Vector{Float64}(undef, N)

sig_mat_abs   = Vector{Float64}(undef, N)
sig_mat_phase = Vector{Float64}(undef, N)
sig_mat_Mz    = Vector{Float64}(undef, N)

# Simulation
for (i, spin) ∈ enumerate(spins)
    ss_vec_Mz[i] = spin.M_ss[3]
    signal_vec = complex(spin.M_ss[1], spin.M_ss[2])
    ss_vec_abs[i] = abs(signal_vec)
    ss_vec_phase[i] = angle(signal_vec)

    ss = steady_state(spin)
    sig_abs[i]   = abs(ss)
    sig_phase[i] = angle(ss)
    
    ss_geo = steady_state_geometric(spin)
    sig_geo_abs[i] = abs(ss_geo)
    sig_geo_Mz[i]  = steady_state_geometric_Mz(spin)

    
    ss_mat = steady_state_matrix(spin)
    signal = complex(ss_mat.x, ss_mat.y)
    sig_mat_abs[i] = abs(signal)
    sig_mat_phase[i] = angle(signal)
    sig_mat_Mz[i] = ss_mat.z
end

# Plotting results
pMag   = plot(B0_vals, sig_mat_abs, label = false, lw = 2.5, color = 4,
         xlabel = "Offset [Hz]", ylabel = "Magnitude", title = "bSSFP Off-Resonace profile")
         #plot!(pMag, B0_vals, sig_geo_abs, label = "Geometric", lw= 2)
         #plot!(pMag, B0_vals, sig_abs, label = "Manual", lw= 2) 
         #plot!(pMag, B0_vals, ss_vec_abs, label = "Spin calculated", lw= 2)
pPhase = plot(B0_vals, sig_mat_phase, label = false, ylabel = "Phase [rad]", lw= 2.5)
         #plot!(pPhase, B0_vals, sig_phase, label = "Manual", lw= 2)
         #plot!(pPhase, B0_vals, ss_vec_phase, label = "Spin calculated", lw= 2)
pSig   = plot(pMag, pPhase; layout = (2,1))

pMz = plot(B0_vals, sig_mat_Mz, label = false, lw= 2.5,
        xlabel = "Offset [Hz]", ylabel = "Mz - Magnetization", title = "bSSFP Off-Resonace profile")
       #plot!(pMz, B0_vals, sig_geo_Mz, label = "Geometric", lw= 2)
       #plot!(pMz, B0_vals, ss_vec_Mz, label = "Spin calculated", lw= 2)

display(pSig)
# display(pMz)

