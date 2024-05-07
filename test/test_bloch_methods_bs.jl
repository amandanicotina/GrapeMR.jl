using Plots
using GrapeMR
using BlochSim

# """Relaxation: ODE, BlochSim, GrapeMR"""
# # RF Construction
# N = 500;
# rf_time = 500; #[ms]
# Δt = rf_time/N;
# ΔB1, Bz = 0.0, zeros(1,N);
# B1x = zeros(1,N);
# B1y = zeros(1,N);
# # Spin
# Mx₀, My₀, Mz₀ = 0.0, 1.0, 0.0;
# T1, T2, Δf    = 500.0, 300.0, 10.0;
# target = ["-"];
# label  = ["-"];

# (spins, field_init) = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], [Δf], target, label, rf_time*1e-3, B1x, B1y, ΔB1, Bz);
# spin_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf)
# spin_grape = spins[1]

# M = Vector{BlochSim.Magnetization{Float64}}(undef, N + 1);
# M[1] = copy(spin_bs.M);

# # GrapeMR Simulation
#     mag_grape = GrapeMR.forward_propagation(field_init, spin_grape);

# # BlochSim Simulation
#     for n ∈ 1:N
#         (A_fp, B_fp) = BlochSim.freeprecess(spin_bs, Δt)
#         M[n+1] = A_fp * M[n] + B_fp
#     end
#     Mx_bs = getproperty.(M, :x); My_bs = getproperty.(M, :y); Mz_bs = getproperty.(M, :z);
   
#     spin_bs_end = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf);
#     BlochSim.freeprecess!(spin_bs_end, rf_time);

# # ODE Solution
#     time = range(0.0, (rf_time*1e-3), N+1);
#     θ = 2π*(time.*spin_grape.B0inho);
#     Mx_sol = (Mx₀*cos.(θ) + My₀*sin.(θ)).*exp.(-spin_grape.T2.*time);
#     My_sol = (My₀*cos.(θ) - Mx₀*sin.(θ)).*exp.(-spin_grape.T2.*time);
#     Mz_sol = 1.0 .- (1.0 .- spin_grape.M_init[3])*exp.(-spin_grape.T1.*time);

# # Plots
#     time = range(0.0, rf_time, N+1)
#     pMx = plot(time, Mx_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mx")
#         plot!(pMx, time, mag_grape[2,:], label = "Forward Propagation")
#         plot!(pMx, time, Mx_sol, label = "ODE Solution")
#         #scatter!([time[end]], [spin_bs_end.M.x], label = "BlochSim final val")

#     pMy = plot(time, My_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
#         plot!(pMy, time, mag_grape[3,:], label = "Forward Propagation")
#         plot!(pMy, time, My_sol, label = "ODE Solution")
#         #scatter!([time[end]], [spin_bs_end.M.y], label = "BlochSim final val")

#     pMz = plot(time, Mz_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
#         plot!(pMz, time, mag_grape[4,:], label = "Forward Propagation")
#         plot!(pMz, time, Mz_sol, label = "ODE Solution")
#         #scatter!([time[end]], [spin_bs_end.M.z], label = "BlochSim final val")
#     display(pMz)
#     pMag = plot(pMx, pMy; layout= (2,1))

"""Rotation: ODE, BlochSim, GrapeMR"""
# RF Construction
    N = 500;
    rf_time = 1000; #[ms]
    Δt = rf_time/N;
    ΔB1, Bz = 0.0, zeros(1,N);
    B1x = 5*ones(1,N);
    B1y = zeros(1,N);
# Spin
    Mx₀, My₀, Mz₀ = 0.0, 0.0, 1.0;
    T1, T2, Δf    = 500000.0, 300000.0, 0.0;
    target = ["-"];
    label  = ["-"];

    (spins, field_init) = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], [Δf], target, label, rf_time*1e-3, B1x, B1y, ΔB1, Bz);
    spin_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf)
    spin_grape = spins[1]

    M = Vector{BlochSim.Magnetization{Float64}}(undef, N + 1);
    M[1] = copy(spin_bs.M);

# GrapeMR Simulation
    mag_grape = GrapeMR.forward_propagation(field_init, spin_grape);

# BlochSim Simulation
    waveform_G = vec(B1x).*(1e4/γ_¹H);
    rf_full = BlochSim.RF(waveform_G, Δt);
    rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G];

    pRF = plot(round.(vec(B1x).*(2π*Δt*1e-3), digits = 4), label = "GrapeMR")
            scatter!(pRF, round.(rf_full.α, digits = 4), label = "BlochSim")
    for n ∈ 1:N
        (A, B) = BlochSim.excite(spin_bs, rf_disc[n])
        M[n+1] = A * M[n] + B
    end
    Mx_bs = getproperty.(M, :x); My_bs = getproperty.(M, :y); Mz_bs = getproperty.(M, :z);

# Plots
    time = range(0.0, rf_time, N+1)
    pMx = plot(time, Mx_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mx")
        plot!(pMx, time, mag_grape[2,:], label = "Forward Propagation")

    pMy = plot(time, My_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
        plot!(pMy, time, mag_grape[3,:], label = "Forward Propagation")

    pMz = plot(time, Mz_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
        plot!(pMz, time, mag_grape[4,:], label = "Forward Propagation")