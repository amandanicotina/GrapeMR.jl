using Plots
using GrapeMR
using BlochSim
using ParameterSchedulers

α_rf_grape = [];
function test_rf(wave_Hz, Δt)
    α_rf_grape = wave_Hz.*(2π*Δt*1e-3)

end


# Spin
    M0_bs, T1, T2, Δf = 1.0, 500, 300, 0; #[ms] and [Hz]
    M0 = [0.0; 0.0; 1.0];
    target = ["-"];
    label = ["-"];

# RF Construction
    N = 500;
    rf_time = 1000; #[ms]
    Δt = rf_time/N; #[ms]

    B0, ΔB1, Bz = [0.0], [1.0], zeros(1,N);
    B1x = 10*ones(1,N);
    B1y = zeros(1,N); 

    ##### NORMALIZE #####
    (spins, field_init) = normalization(M0, [T1*1e-3], [T2*1e-3], B0, target, label, rf_time*1e-3, B1x, B1y, ΔB1, Bz);

    ##### OPTIMIZE #####
    max_iter = 300
    lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1) 
    opt_params   = OptimizationParams(N, max_iter, cost_target_one_spin, [false false false]);
    grape_output = GrapeMR.grape(opt_params, field_init, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 

    waveform_Hz = vec(grape_output.control_field.B1x); #5*ones(N);
    waveform_G  = waveform_Hz.*(1e4/γ_¹H);

# BlochSim Simulation #
    rf_full = BlochSim.RF(waveform_G, Δt)
    rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G]
    
    pRF = plot(waveform_Hz.*(2π*Δt*1e-3), label = "GrapeMR")
            plot!(pRF, rf_full.α, label = "BlochSim")

    # Spin Construction
    spin = BlochSim.Spin(M0_bs, T1, T2, Δf);
    M = Vector{BlochSim.Magnetization{Float64}}(undef, N + 1);
    M[1] = copy(spin.M);

    # Simulation
    for n ∈ 1:N
        (A, B) = BlochSim.excite(spin, rf_disc[n])
        M[n+1] = A * M[n] + B
    end
    Mx = getproperty.(M, :x); My = getproperty.(M, :y); Mz = getproperty.(M, :z);

# GrapeMR Simulation #
    grape_mag_opt = grape_output.isochromats[1].magnetization.dynamics;
    mag_test  = forward_propagation(field_init, spins[1])

# Plots
    time = range(0.0, rf_time, N+1)
    pMx = plot(time, My, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mx")
        plot!(pMx, time, mag_test[2,:], label = "Forward Propagation")
        # scatter!(pMx, time, grape_mag_opt[2,:], label = "GrapeMR")

    pMy = plot(time, Mx, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
        plot!(pMy, time, mag_test[3,:], label = "Forward Propagation")
        # scatter!(pMy, time, grape_mag_opt[3,:], label = "GrapeMR")

    pMz = plot(time, Mz, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
        plot!(pMz, time, mag_test[4,:], label = "Forward Propagation")
       # scatter!(pMz, time, grape_mag_opt[4,:], label = "GrapeMR")
