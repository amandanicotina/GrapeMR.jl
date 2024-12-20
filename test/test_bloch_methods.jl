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

function test_relaxation_dynamics()
    # Parameters
    N = 1000;
    rf_time = 1000.0; #[ms]
    Δt = rf_time/N;
    # Spin
    Mx₀, My₀, Mz₀ = 0.0, 1.0, 0.0;
    T1, T2, Δf    = 500.0, 300.0, [0.0, -5.0, -30.0];
    target = ["-", "-", "-"];
    label  = [string("Δf = ", Δf[1]), string("Δf = ", Δf[2]), string("Δf = ", Δf[3])];
    # RFs
    ΔB1, Bz = 1.0, zeros(1,N), 0.0;
    Bx_off = zeros(1,N); 
    By_off = zeros(1,N); 
    B_ref  = 1π

    # Control Field and Spins objects
    control_field = ControlField(Bx_off, By_off, B_ref, Bz, rf_time*1e-3)
    spins         = GrapeMR.Spin([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], Δf, ΔB1, target, label)

    # Constructing different Spin objects to be tested
    spins_bs = BlochSim.Spin.(Ref(BlochSim.Magnetization(Mx₀, My₀, Mz₀)), 1, T1, T2, Δf)

    # Magnetization dynamics for all spin objects
    mag_grape = GrapeMR.forward_propagation.(Ref(control_field), spins)
    mag_bs    = BlochSim_simulation.(spins_bs, Δt, N)

    return mag_bs, mag_grape, rf_time
end

function test_rotation_dynamics()
    # Parameters
    N = 1000;
    rf_time = 500; #[ms]
    Δt = rf_time/N
    # Spin
    Mx₀, My₀, Mz₀ = 0.0, 0.0, 1.0;
    T1, T2,       = [1e10], [1e10];
    Δf = [0.0, -5.0, -30.0];
    target = ["[0.0, 1.0, 0.0]", "[0.0, 1.0, 0.0]", "[0.0, 1.0, 0.0]"];
    label  = [string("Δf = ", Δf[1]), string("Δf = ", Δf[2]), string("Δf = ", Δf[3])];
    # RFs
    ΔB1, Bz = 1.0, zeros(1,N), 0.0;
    B1ref = 1.0

    # Control Field and Spins objects
    control_field = hard_RF(N, rf_time*1e-3, B1ref)
    spins         = GrapeMR.Spin([Mx₀, My₀, Mz₀], T1, T2, Δf, ΔB1, target, label)

    # Constructing different Spin objects to be tested ; 
    spins_bs = BlochSim.Spin.(Ref(BlochSim.Magnetization(Mx₀, My₀, Mz₀)), 1, T1, T2, Δf)

    # Convert RF for BlochSim
    waveform_T = control_field.B1x./γ_¹H .+ im*control_field.B1y./γ_¹H;
    waveform_G = vec(waveform_T).*1e4;
    rf_full = BlochSim.RF(waveform_G, Δt)
    rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G]

    # Calling functions for all objects
    mag_grape = GrapeMR.forward_propagation.(Ref(control_field), spins);
    mag_bs    = BlochSim_simulation.(spins_bs, Ref(rf_disc), N)

    return mag_bs, mag_grape, rf_time
end

function test_shaped_pulse_dynamics()
    # Parameters
    N = 1000;
    rf_time = 1000; #[ms]
    Δt = rf_time/N;

    # Spin
    Mx₀, My₀, Mz₀ = 0.0, 0.0, 1.0;
    T1, T2, Δf, ΔB1 = 500.0, 300.0, [0.0, -5.0, -30.0], [1.0];
    target = ["[0.0, 1.0, 0.0]"];
    label  = [string("Δf = ", Δf[1]), string("Δf = ", Δf[2]), string("Δf = ", Δf[3])];
    spins  = GrapeMR.Spin([Mx₀, My₀, Mz₀], [T1*1e-3], [T2*1e-3], Δf, ΔB1, target, label)

    # BlochSim spin objects
    spins_bs = BlochSim.Spin.(Ref(BlochSim.Magnetization(Mx₀, My₀, Mz₀)), 1, T1, T2, Δf)

    # Parameters
    max_iter     = 1000; start=1e-1; degree=2; max_iter=max_iter+1
    opt_params   = OptimizationParams(start, degree, Int(ceil(max_iter)))
    grape_params = GrapeParams(N, GrapeMR.euclidean_norm, Dict("B1x" => true, "B1y" => true, "Bz" => false))
    params       = Parameters(grape_params, opt_params)

    # Initial RF Pulse
    B1ref = 1.0
    control_field = spline_RF(grape_params.N, rf_time*1e-3, B1ref)

    # Run Optimization
    grape_output = @time grape(params, control_field, spins); 

    # BlochSim Simulation 
    waveform_Hz = grape_output.control_field.B1x .- im*grape_output.control_field.B1y;
    waveform_T = grape_output.control_field.B1x./γ_¹H .- im*grape_output.control_field.B1y./γ_¹H;
    waveform_G = vec(waveform_T).*1e4;
    rf_full = BlochSim.RF(waveform_G, Δt)
    rf_disc = [BlochSim.RF([Δrf], Δt) for Δrf ∈ waveform_G]

    pRF = plot(abs.(waveform_Hz)'.*(2π*Δt*1e-3), label = "GrapeMR")
            plot!(pRF, rf_full.α, label = "BlochSim")

    # Spin Construction
    mag_bs_opt   = BlochSim_simulation.(spins_bs, Ref(rf_disc), N)
    mag_gp_opt   = [iso.magnetization.dynamics for iso ∈ grape_output.isochromats]

    return mag_bs_opt, mag_gp_opt, rf_time
end

function plot_magnetization_tracking(M_bs::Vector{BlochSim.Magnetization{Float64}}, M_gp::Matrix{Float64}, t_rf::Real)
    Mx_bs, My_bs, Mz_bs    = getproperty.(M_bs, :x), getproperty.(M_bs, :y), getproperty.(M_bs, :z);
    Mx_gp, My_gp, Mz_gp    = M_gp[2,:], M_gp[3,:], M_gp[4,:]

    time = range(0.0, stop=t_rf, length=length(Mx_gp))

    pMx = plot(time, -My_bs, label = "BlochSim", ylabel = "Magnetization", title = "Mx")
        plot!(pMx, time, Mx_gp, label = "Forward Propagation")

    pMy = plot(time, Mx_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "My")
        plot!(pMy, time, My_gp, label = "Forward Propagation")

    pMz = plot(time, Mz_bs, label = "BlochSim", xlabel = "t (ms)", ylabel = "Magnetization", title = "Mz")
        plot!(pMz, time, Mz_gp, label = "Forward Propagation")
    
    pMag = plot(pMx, pMy; layout= (2,1), xlims = [0.0, t_rf], ylims = [-1.05, 1.05])

    display(pMag)
    display(pMz)
    return pMag, pMz
end 

# -----------------------------------------------------------------------
# Test Rotation Dynamics
# -----------------------------------------------------------------------
(mag_bs_rot, mag_gp_rot, rf_time_rot)  = test_rotation_dynamics()
Mx_bs, My_bs, Mz_bs    = getproperty.(mag_bs_rot[1], :x), getproperty.(mag_bs_rot[1], :y), getproperty.(mag_bs_rot[1], :z);
Mx_gp, My_gp, Mz_gp    = mag_gp_rot[1][2,:], mag_gp_rot[1][3,:], mag_gp_rot[1][4,:]

@test all(round.(Mx_gp, digits = 1) .== round.(My_bs, digits = 1))
@test all(round.(My_gp, digits = 1) .== round.(Mx_bs, digits = 1))
@test all(round.(Mz_gp, digits = 1) .== round.(Mz_bs, digits = 1))

# -----------------------------------------------------------------------
# Test Relaxation Dynamics
# -----------------------------------------------------------------------
(mag_bs_relax, mag_gp_relax, rf_time_relax) = test_relaxation_dynamics()
Mx_bs_relax, My_bs_relax, Mz_bs_relax = getproperty.(mag_bs_relax[1], :x), getproperty.(mag_bs_relax[1], :y), getproperty.(mag_bs_relax[1], :z);
Mx_gp_relax, My_gp_relax, Mz_gp_relax = mag_gp_relax[1][2,:], mag_gp_relax[1][3,:], mag_gp_relax[1][4,:]

@test all(round.(Mx_gp_relax, digits = 1) .== round.(Mx_bs_relax, digits = 1))
@test all(round.(My_gp_relax, digits = 1) .== round.(My_bs_relax, digits = 1))
@test all(round.(Mz_gp_relax, digits = 1) .== round.(Mz_bs_relax, digits = 1))

# -----------------------------------------------------------------------
# Test Shaped Pulse Dynamics
# -----------------------------------------------------------------------
(mag_bs_shape, mag_gp_shape, rf_time_shape) = test_shaped_pulse_dynamics()
Mx_bs_shape, My_bs_shape, Mz_bs_shape = getproperty.(mag_bs_shape[1], :x), getproperty.(mag_bs_shape[1], :y), getproperty.(mag_bs_shape[1], :z);
Mx_gp_shape, My_gp_shape, Mz_gp_shape = mag_gp_shape[1][2,:], mag_gp_shape[1][3,:], mag_gp_shape[1][4,:]

# tol = 1e-1  
# @test all(isapprox.(Mx_gp_shape, Mx_bs_shape, atol = tol))
# @test all(isapprox.(My_gp_shape, My_bs_shape, atol = tol))
# @test all(isapprox.(Mz_gp_shape, Mz_bs_shape, atol = tol))  

# @test all(round.(Mx_gp_shape, digits = 3) .== round.(My_bs_shape, digits = 3))
# @test all(round.(My_gp_shape, digits = 3) .== round.(Mx_bs_shape, digits = 3))
# @test all(round.(Mz_gp_shape, digits = 3) .== round.(Mz_bs_shape, digits = 3))