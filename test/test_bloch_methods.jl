using Plots
using GrapeMR
using BlochSim
using ParameterSchedulers

# # TODO actually test each component of the magnetization not just ploting it -> make it into a real test
# # TODO figure out a better way to check for all different combinations

# ODE Solution
function ODE_simulation(spin::GrapeMR.Spin, t_rf::Real, N::Int)
    T₁, T₂ = spin.T1[], spin.T2[]
    Mx₀, My₀, Mz₀ = spin.M_init[1], spin.M_init[2], spin.M_init[3]

    time = range(0.0, stop=(t_rf*1e-3), length=N+1)
    θ = 2π * (time * spin.B0inho)
    cosθ = cos.(θ)
    sinθ = sin.(θ)
    expT2 = exp.(-T₂ * time)
    expT1 = exp.(-T₁ * time)

    M_sol = zeros(3, length(time))  # Preallocate a 3xN matrix
    M_sol[1, :] = (Mx₀ * cosθ + My₀ * sinθ) .* expT2
    M_sol[2, :] = (My₀ * cosθ - Mx₀ * sinθ) .* expT2
    M_sol[3, :] = 1.0 .- (1.0 - Mz₀) .* expT1

    return M_sol
end

# function ODE_simulation(spin::GrapeMR.Spin, rf, N::Int)
#     Mx₀, My₀, Mz₀ = spin.M_init[1], spin.M_init[2], spin.M_init[3]
#     M0 = [Mx₀, My₀, Mz₀];
#     M_sol = zeros(3, length(time)) 
#     M_sol[:,1] = M0
    
#     B1 = sqrt((rf.B1x.^2) + (rf.B1y.^2))
#     time = range(0.0, stop=(t_rf*1e-3), length=N+1)
#     α_arr = 2π * (time .* B1)

#     for (i, α) ∈ enumerate(α_arr)
#         cosα = cos.(α)
#         sinα = sin.(α)

#         # Rotation matrices
#         Rx = [1.0  0.0     0.0;
#             0.0  cosα  -sinα;
#             0.0  sinα  cosα]

#         Ry = [cosα  0.0   sinα;
#             0.0     1.0   0.0;
#             -sinα 0.0   cosα]
        
#         M_sol[:, i+1] = Rx*M[:,i]
#     end

#     return M_sol
# end


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
function plot_magnetization_tracking(M_bs::Vector{BlochSim.Magnetization{Float64}}, M_gp::Matrix{Float64}, M_ODE::Matrix{Float64}, t_rf::Real)
    Mx_bs, My_bs, Mz_bs    = getproperty.(M_bs, :x), getproperty.(M_bs, :y), getproperty.(M_bs, :z);
    Mx_gp, My_gp, Mz_gp    = M_gp[2,:], M_gp[3,:], M_gp[4,:]
    Mx_ODE, My_ODE, Mz_ODE = M_ODE[1,:], M_ODE[2,:], M_ODE[3,:]

    time = range(0.0, stop=t_rf, length=length(Mx_gp))

    pMx = plot(time, Mx_bs, label = "BlochSim", ylabel = "Magnetization", title = "Mx")
        plot!(pMx, time, Mx_gp, label = "Forward Propagation")
        #plot!(pMx, time, Mx_ODE, label = "ODE Solution")
        #scatter!([time[end]], [spin_bs_end.M.x], label = "BlochSim final val")

    pMy = plot(time, My_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
        plot!(pMy, time, My_gp, label = "Forward Propagation")
        #plot!(pMy, time, My_ODE, label = "ODE Solution")
        #scatter!([time[end]], [spin_bs_end.M.y], label = "BlochSim final val")

    pMz = plot(time, Mz_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
        plot!(pMz, time, Mz_gp, label = "Forward Propagation")
        #plot!(pMz, time, Mz_ODE, label = "ODE Solution")
        #scatter!([time[end]], [spin_bs_end.M.z], label = "BlochSim final val")
    
    pMag = plot(pMx, pMy; layout= (2,1), xlims = [0.0, t_rf], ylims = [-1.05, 1.05])
    return pMag, pMz
end

###################### TESTS ######################

# @safetestset begin "Free Precession/Relaxation: ODE vs BlochSim vs Forward Propagation"
# Parameters
# N = 500;
# rf_time = 1000; #[ms]
# Δt = rf_time/N;
# # Spin
# Mx₀, My₀, Mz₀ = 0.0, 1.0, 0.0;
# T1, T2,       = 500.0, 300.0;
# Δf = [0.0, -5.0, -30.0];
# target = ["-"];
# label  = ["-"];
# # RFs
# ΔB1, Bz = 1.0, zeros(1,N), 0.0;
# Bx_off  = zeros(1,N); 
# By_off  = zeros(1,N); 

# # Constructing different Spin objects to be tested
# (spins, control_field) = normalization([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], Δf, target, label, rf_time*1e-3, Bx_off, By_off, ΔB1, Bz); 
# spin1_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[1])
# spin2_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[2])
# spin3_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[3])

# # Calling functions for all objects
# #for (i, B0) ∈ enumerate(Δf)
# mag_grape = GrapeMR.forward_propagation(control_field, spins[3])
# mag_bs    = BlochSim_simulation(spin3_bs, Δt, N)
# mag_ODE   = ODE_simulation(spins[3], rf_time, N)
# #end

# (pMag, pMz) = plot_magnetization_tracking(mag_bs, mag_grape, mag_ODE, rf_time)
# display(pMz)
# display(pMag) 
# end

"""Rotation: ODE vs BlochSim vs Forward Propagation"""
# @safetestset begin "Rotation: ODE vs BlochSim vs Forward Propagation"
# Parameters
N = 500;
rf_time = 500; #[ms]
Δt = rf_time/N;
# Spin
Mx₀, My₀, Mz₀ = 0.0, 0.0, 1.0;
T1 = 500.0;
T2 = 300.0;
Δf = [0.0, -5.0, 8.0];
target = ["-"];
label  = ["-"];
# RFs
ΔB1, Bz = 1.0, zeros(1,N), 0.0;
B1x = 8*initial_field_spline(N, rf_time*1e-3)'; # rand(1, N); #
B1y = zeros(1,N); #5*initial_field_spline(N, rf_time*1e-3)'; # 
control_field = ControlField(B1x, B1y, 1.0, Bz, rf_time*1e-3)

# Constructing different Spin objects to be tested 
sp = GrapeMR.Spin([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], Δf, ΔB1, target, label)
spin1_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[1])
spin2_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[2])
spin3_bs = BlochSim.Spin(BlochSim.Magnetization(Mx₀, My₀, Mz₀), 1, T1, T2, Δf[3])

max_iter = 300
lr_scheduler = Poly(start=1e-1, degree=2, max_iter=max_iter+1) 
opt_params   = OptimizationParams(N, max_iter, cost_target_one_spin, [false false false]);
grape_output = GrapeMR.grape(opt_params, control_field, spins, lr_scheduler; max_iter = max_iter, ϵ=1e-2); 

# Convert RF for BlochSim
waveform_T = grape_output.control_field.B1x./γ_¹H .+ im*grape_output.control_field.B1y./γ_¹H;
waveform_G = vec(waveform_T).*1e4;
rf_full = BlochSim.RF(waveform_G, Δt)
rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G]

# Calling functions for all objects
mag_grape = GrapeMR.forward_propagation(control_field, sp[3]);
mag_bs    = BlochSim_simulation(spin3_bs, rf_disc, N)
mag_ODE   = ODE_simulation(spins[1], rf_time, N);

(pMag, pMz) = plot_magnetization_tracking(mag_bs, mag_grape, mag_ODE, rf_time)
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