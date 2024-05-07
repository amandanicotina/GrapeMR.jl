using Plots
using GrapeMR
using BlochSim
using ParameterSchedulers

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

"""Rotation: BlochSim, GrapeMR in units, GrapeMR normalized"""
# RF Construction
    N = 500;
    rf_time = 500; #[ms]
    Δt = rf_time/N;
    ΔB1, Bz = 1.0, zeros(1,N);
    B1x = 3*initial_field_spline(N, rf_time*1e-3)'; 
    B1y = 5*initial_field_spline(N, rf_time*1e-3)'; 
    B_ref = maximum(B1x);

# Spin
    Mx₀, My₀, Mz₀ = 0.0, 0.0, 1.0;
    T1, T2, Δf    = 500.0, 300.0, 0.0;
    target = ["-"];
    label  = ["-"];

# GrapeMR Simulation
    # Normalized    
    (spins, field_init) = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], [Δf], target, label, rf_time*1e-3, B1x, B1y, ΔB1, Bz);
    spin_grape = spins[1];
    
    mag_grape_norm = GrapeMR.forward_propagation(grape_output.control_field, spin_grape);

    # Rescale with GrapeMR
    max_iter = 300;
    lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1); 
    opt_params   = OptimizationParams(N, max_iter, cost_target_one_spin, [true true false]);
    grape_output = GrapeMR.grape(opt_params, field_init, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 
   
    mag_grape_units = grape_output.isochromats[1].magnetization.dynamics;
    mag_grape_norm = GrapeMR.forward_propagation(grape_output.control_field, spin_grape);

# BlochSim Simulation
    spin_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf)
    M = Vector{BlochSim.Magnetization{Float64}}(undef, N + 1);
    M[1] = copy(spin_bs.M);

    waveform_T = grape_output.control_field.B1x./γ_¹H .+ im*grape_output.control_field.B1y./γ_¹H
    waveform_G = vec(waveform_T).*1e4;
    rf_full = BlochSim.RF(waveform_G, Δt);
    rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G];

    for n ∈ 1:N
        (A, B) = BlochSim.excite(spin_bs, rf_disc[n])
        M[n+1] = A * M[n] + B
    end
    Mx_bs = round.(getproperty.(M, :x), digits = 5); My_bs = round.(getproperty.(M, :y), digits = 5); Mz_bs = round.(getproperty.(M, :z), digits = 5);

# Plots
    waveform = vec(γ_¹H.*waveform_T.*(2π*Δt*1e-3))
    pFlip = plot(round.(real.(waveform), digits = 4), label = "GrapeMR")
        plot!(pFlip, round.(rf_full.α, digits = 4), label = "BlochSim")

    pPhase = plot(round.(angle.(waveform), digits = 4), label = "GrapeMR")
        plot!(pPhase, round.(rf_full.θ, digits = 4), label = "BlochSim")
    
    pRF = plot(pFlip, pPhase; layout= (2,1))
    
    time = range(0.0, rf_time, N+1)
    pMx = plot(time, My_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mx")
        plot!(pMx, time, mag_grape_norm[2,:], label = "Forward Propagation")
        plot!(pMx, time, mag_grape_units[2,:], label = "GrapeMR Units")

    pMy = plot(time, Mx_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
        plot!(pMy, time, mag_grape_norm[3,:], label = "Forward Propagation")
        plot!(pMy, time, mag_grape_units[3,:], label = "GrapeMR Units")

    pMz = plot(time, Mz_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
        plot!(pMz, time, mag_grape_norm[4,:], label = "Forward Propagation")
        plot!(pMz, time, mag_grape_units[4,:], label = "GrapeMR Units")
        
    pMag = plot(pMx, pMy, pMz; layout= (3,1))