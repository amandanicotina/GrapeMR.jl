using Plots
using GrapeMR
using BlochSim

# BlochSim Simulation: Free Precession/Relaxation
function BlochSim_simulation(spin::BlochSim.Spin, Δt::Float64, N::Int)
    M = Vector{BlochSim.Magnetization{Float64}}(undef, N + 1);
    M[1] = copy(spin.M);
    for n ∈ 1:N
        (A_fp, B_fp) = BlochSim.freeprecess(spin, Δt)
        M[n+1] = A_fp * M[n] + B_fp
    end
    return M
end

# BlochSim Simulation: Rotation
function BlochSim_simulation(spin::BlochSim.Spin, rf::Vector{RF{Float64, Gradient{Int64}}}, N::Int)
    M = Vector{BlochSim.Magnetization{Float64}}(undef, N + 1);
    M[1] = copy(spin.M);
    for n ∈ 1:N
        (A_ex, B_ex) = BlochSim.excite(spin, rf[n])
        M[n+1] = A_ex * M[n] + B_ex
    end
    return M
end

# Plots
function plot_magnetization_tracking(M_bs::Vector{BlochSim.Magnetization{Float64}}, M_gp::Matrix{Float64}, t_rf::Real)
    Mx_bs, My_bs, Mz_bs    = getproperty.(M_bs, :x), getproperty.(M_bs, :y), getproperty.(M_bs, :z);
    Mx_gp, My_gp, Mz_gp    = M_gp[2,:], M_gp[3,:], M_gp[4,:]

    time = range(0.0, stop=t_rf, length=length(Mx_gp))

    pMx = plot(time, Mx_bs, label = "BlochSim", ylabel = "Magnetization", title = "Mx")
        plot!(pMx, time, Mx_gp, label = "Forward Propagation")


    pMy = plot(time, My_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
        plot!(pMy, time, My_gp, label = "Forward Propagation")

    pMz = plot(time, Mz_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
        plot!(pMz, time, Mz_gp, label = "Forward Propagation")
    
    pMag = plot(pMx, pMy; layout= (2,1), xlims = [0.0, t_rf], ylims = [-1.05, 1.05])
    return pMag, pMz
end

###################### TESTS ######################
# @safetestset begin "Free Precession/Relaxation: ODE vs BlochSim vs Forward Propagation"
# Parameters
N = 1000;
rf_time = 1000; #[ms]
Δt = rf_time/N;
# Spin
Mx₀, My₀, Mz₀ = 0.0, 1.0, 0.0;
T1, T2,       = 500.0, 300.0;
Δf = [0.0, -5.0, -30.0];
target = ["-", "-", "-"];
label  = [string("Δf = ", Δf[1]), string("Δf = ", Δf[2]), string("Δf = ", Δf[3])];
# RFs
ΔB1, Bz = 1.0, zeros(1,N), 0.0;
Bx_off  = zeros(1,N); 
By_off  = zeros(1,N); 
B_ref = π

# Control Field and Spins objects
control_field = ControlField(Bx_off, By_off, B_ref, Bz, rf_time*1e-3)
spins         = GrapeMR.Spin([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], Δf, ΔB1, target, label)

# Constructing different Spin objects to be tested
(norm_spins, norm_control_field) = normalization(spins, control_field); 
spin1_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[1])
spin2_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[2])
spin3_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[3])

# Calling functions for all objects
#for (i, B0) ∈ enumerate(Δf)
mag_grape = GrapeMR.forward_propagation(norm_control_field, norm_spins[1])
mag_bs    = BlochSim_simulation(spin1_bs, Δt, N)
#end

(pMag, pMz) = plot_magnetization_tracking(mag_bs, mag_grape, rf_time)
display(pMz)
display(pMag) 
# end

"""Rotation: ODE vs BlochSim vs Forward Propagation"""
# @safetestset begin "Rotation: ODE vs BlochSim vs Forward Propagation"
# Parameters
N = 500;
rf_time = 500; #[ms]
Δt = rf_time/N;
# Spin
Mx₀, My₀, Mz₀ = 0.0, 1.0, 0.0;
T1, T2,       = 500.0, 300.0;
Δf = [0.0, -5.0, -30.0];
target = ["[0.0, 1.0, 0.0]", "[0.0, 1.0, 0.0]", "[0.0, 1.0, 0.0]"];
label  = [string("Δf = ", Δf[1]), string("Δf = ", Δf[2]), string("Δf = ", Δf[3])];
# RFs
ΔB1, Bz = 1.0, zeros(1,N), 0.0;
B1_ref = 1.0
B1x = spline_RF(N, rf_time*1e-3)'; # rand(1, N); #
B1y = zeros(1,N); #5*initial_field_spline(N, rf_time*1e-3)'; # 

# Control Field and Spins objects
control_field = ControlField(B1x, B1y, B1_ref, Bz, rf_time*1e-3)
spins         = GrapeMR.Spin([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], Δf, ΔB1, target, label)

# Constructing different Spin objects to be tested 
(norm_spins, norm_control_field) = normalization(spins, control_field); 
spin1_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[1])
spin2_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[2])
spin3_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[3])

max_iter = 300
start=1e-1; degree=2; max_iter=max_iter+1
opt_params   = OptimizationParams(start, degree, Int(ceil(max_iter)))
grape_params = GrapeParams(N, :target_one_spin, [false false false])
params = Parameters(grape_params, opt_params)
grape_output = norm_grape(params, control_field, spins);
(back_spins, back_cf) = inverse_normalization(grape_output.isochromats.spins, grape_output.control_field)


# Convert RF for BlochSim
waveform_T = norm_control_field.B1x .+ im*norm_control_field.B1y./γ_¹H;
waveform_G = vec(waveform_T).*1e4;
rf_full = BlochSim.RF(waveform_G, Δt)
rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G]

# Calling functions for all objects
mag_grape = GrapeMR.forward_propagation(norm_control_field, norm_spins[1]);
mag_bs    = BlochSim_simulation(spin1_bs, rf_disc, N)

(pMag, pMz) = plot_magnetization_tracking(mag_bs, mag_grape, rf_time)
display(pMag)
display(pMz) 



# function RF_testing(rf::ControlField, Δt::Float64)
#     waveform_T = rf.B1x./γ_¹H .+ im*rf.B1y./γ_¹H
#     waveform_G = vec(waveform_T).*1e4;

#     rf_BS = BlochSim.RF(waveform_G, Δt);
#     rf_grape = vec(γ_¹H.*waveform_T.*(2π*Δt*1e-3));
    
#     steps = length(rf.B1x)
#     t = range(0.0, Δt*steps, steps);

#     return rf_BS, rf_grape, t
# end
# rf_data = RF_testing(control_field, Δt)
# function plot_RF_testing(rf_data)
#     rf_data_BS    = rf_data[1];
#     rf_data_grape = rf_data[2];
#     rf_data_time  = rf_data[3];

#     pFlip = plot(rf_data_time, round.(abs.(rf_data_grape), digits = 4), label = "GrapeMR", xlabel = "t[ms]", ylabel = "Flip angle [rads]", title = "RF pulse")
#         plot!(pFlip, rf_data_time, round.(rf_data_BS.α, digits = 4), label = "BlochSim")

#     pPhase = plot(rf_data_time, round.(angle.(rf_data_grape), digits = 4), label = "GrapeMR", xlabel = "t[ms]", ylabel = "Phase [rads]")
#         plot!(pPhase, rf_data_time, round.(rf_data_BS.θ, digits = 4), label = "BlochSim")
    
#     pRF = plot(pFlip, pPhase; layout= (2,1))

#     return pRF 
# end
# plot_RF_testing(rf_data)


using Plots
using GrapeMR
using BlochSim
using ParameterSchedulers


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
