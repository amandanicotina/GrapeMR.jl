using Plots
using GrapeMR
using BlochSim
using ParameterSchedulers

@safetestset begin "RF pulse: GrapeMR vs BlochSim"
# TODO add optimizated Control Filed ot trest after refactoting
# TODO create script to test each element individually

# Parameters
    N = 500;
    rf_time = 500; #[ms]
    Δt = rf_time/N;
# Spin
    Mx₀, My₀, Mz₀ = 0.0, 0.0, 1.0;
    T1, T2, Δf    = 500.0, 300.0, 0.0;
    target = ["-"];
    label  = ["-"];
# RFs
    # Bx -> on, By -> off
    ΔB1, Bz = 1.0, zeros(1,N), 0.0;
    Bx_on   = 50*initial_field_spline(N, rf_time*1e-3)'; 
    By_off  = zeros(1,N); 
    # Bx -> off, By -> on
    ΔB1, Bz = 1.0, zeros(1,N), 0.0;
    Bx_off  = zeros(1,N);
    By_on   = 7*initial_field_spline(N, rf_time*1e-3)'; 
    # Bx -> on, By -> on
    ΔB1, Bz   = 1.0, zeros(1,N), 0.0;
    Bx_spline = 500*initial_field_spline(N, rf_time*1e-3)'; 
    By_spline = 500*initial_field_spline(N, rf_time*1e-3)'; 

    (-, control_BxOnByOff) = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], [Δf], target, label, rf_time*1e-3, Bx_on, By_off, ΔB1, Bz);
    (-, control_BxOffByOn) = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], [Δf], target, label, rf_time*1e-3, Bx_off, By_on, ΔB1, Bz);
    (-, control_spline)    = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], [Δf], target, label, rf_time*1e-3, Bx_spline, By_spline, ΔB1, Bz);

function RF_testing(rf::ControlField, Δt::Float64)
    waveform_T = rf.B1x./γ_¹H .+ im*rf.B1y./γ_¹H
    waveform_G = vec(waveform_T).*1e4;

    rf_BS = BlochSim.RF(waveform_G, Δt);
    rf_grape = vec(γ_¹H.*waveform_T.*(2π*Δt*1e-3));
    
    steps = length(rf.B1x)
    t = range(0.0, Δt*steps, steps);

    return rf_BS, rf_grape, t
end

# Test RF elements individually
rf_data_BxOnByOff[1]

# Plots
rf_data_BxOnByOff = RF_testing(control_BxOnByOff, Δt)
rf_data_BxOffByOn = RF_testing(control_BxOffByOn, Δt)
rf_data_spline    = RF_testing(control_spline, Δt)

function plot_RF_testing(rf_data)
    rf_data_BS    = rf_data[1];
    rf_data_grape = rf_data[2];
    rf_data_time  = rf_data[3];

    pFlip = plot(rf_data_time, round.(abs.(rf_data_grape), digits = 4), label = "GrapeMR", xlabel = "t[ms]", ylabel = "Flip angle [rads]", title = "RF pulse")
        plot!(pFlip, rf_data_time, round.(rf_data_BS.α, digits = 4), label = "BlochSim")

    pPhase = plot(rf_data_time, round.(angle.(rf_data_grape), digits = 4), label = "GrapeMR", xlabel = "t[ms]", ylabel = "Phase [rads]")
        plot!(pPhase, rf_data_time, round.(rf_data_BS.θ, digits = 4), label = "BlochSim")
    
    pRF = plot(pFlip, pPhase; layout= (2,1))

    return pRF 
end;

    plot_RF_testing(rf_data_BxOnByOff)
    plot_RF_testing(rf_data_BxOffByOn)
    plot_RF_testing(rf_data_spline)
end


