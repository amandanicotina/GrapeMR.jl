using Test, SafeTestsets
using GrapeMR

function test_cost_one_spin(cost_func::Function; target::AbstractVector = [0.0, 1.0, 0.0])
    # Parameters
    N = 1000;
    grape_params = GrapeParams(N, cost_func, Dict("B1x" => true, "B1y" => true, "Bz" => false))

    # RFs
    t_c   = 0.5; #[s]
    B1ref = 1.0

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8];
    T2  = [0.6];
    B0 = 10.0
    offsets = collect(-B0:5:B0)
    ΔB1 = [1.0];
    target = [string(target)];
    label  = ["T1-$(Int(T1[1]*1e3))ms"]

    # Spin and RF objects
    control_field = spline_RF(grape_params.N, t_c, B1ref)
    spins   = GrapeMR.Spin(M0, T1, T2, offsets, ΔB1, target, label)

    iso_vec = dynamics.(control_field, spins)

    # Cost gradient
    cost_vals = grape_params.cost_function.(iso_vec)
    cost_grad = getindex.(cost_vals, 2)
   
    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost.(grape_params.cost_function, iso_vec, ΔM)

    return fd_M, cost_grad
end

function test_cost_steady_state(cost_func::Function)
    # Parameters
    N = 1000;
    grape_params = GrapeParams(N, cost_func, Dict("B1x" => true, "B1y" => true, "Bz" => false))

    # RFs
    t_c   = 0.5; #[s]
    B1ref = 1.0

    # Spin System
    M0 = [0.0, 0.0, 1.0]
    ΔB1 = [1.0]
    B0 = 10.0
    target = ["max"]
    label  = ["s1"]
    T1 = [0.5]
    T2 = [0.1]
    offset = collect(-B0:10:B0)
    α, Δϕ, TR, TE = π/3, 2π, 5e-3, 2.5e-3

    # Spin and RF objects
    control_field = spline_RF(grape_params.N, t_c, B1ref)
    spins = GrapeMR.SteadyState(M0, T1, T2, offset, ΔB1, target, label, α, Δϕ, TR, TE)
    iso_vec = dynamics.(control_field, spins)

    # Cost gradient
    cost_vals = grape_params.cost_function.(iso_vec)
    cost_grad = getindex.(cost_vals, 2)
   
    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost.(grape_params.cost_function, iso_vec, ΔM)

    return fd_M, cost_grad
end

function test_cost_two_spins(cost_func::Function)
    # Parameters
    N = 1000;
    grape_params = GrapeParams(N, cost_func, Dict("B1x" => true, "B1y" => true, "Bz" => false))

    # RFs
    t_c   = 0.5; #[s]
    B1ref = 1.0

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8, 0.4];
    T2  = [0.6, 0.3];
    B0  = [10.0];
    ΔB1 = [1.0];
    target = ["max", "min"];
    label  = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

    # Spin and RF objects
    control_field = spline_RF(grape_params.N, t_c, B1ref)
    spins   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    iso_vec = dynamics.(Ref(control_field), spins)

    # Cost gradient
    cost_vals = grape_params.cost_function.(iso_vec)
    cost_grad = getindex.(cost_vals, 2)

    # Finite difference
    ΔM   = 1e-10;
    fd_M_max = finite_difference_cost(grape_params.cost_function, iso_vec[1], ΔM)
    fd_M_min = finite_difference_cost(grape_params.cost_function, iso_vec[2], ΔM)

    return fd_M_max, fd_M_min, cost_grad
end

(fd_M_en, cost_grad_en) = test_cost_one_spin(GrapeMR.euclidean_norm)
@test all(round.(fd_M_en, digits=3) .== round.(cost_grad_en[:, 2:end], digits=3))

(fd_M_tos, cost_grad_tos) = test_cost_one_spin(GrapeMR.spin_target; target = [1.0, 0.0, 0.0])
@test all(round.(fd_M_tos, digits=3) .== round.(cost_grad_tos[:, 2:end] , digits=3))

(fd_M_ssOffsets, cost_grad_ssOffsets) = test_cost_steady_state(GrapeMR.steady_state_offset_targets)
@test all(round.(fd_M_ssOffsets, digits=3) .== round.(cost_grad_ssOffsets[:, 2:end] , digits=3))

(fd_M_ssSatCon, cost_grad_ssSatCon) = test_cost_steady_state(GrapeMR.saturation_contrast_steady_state)
@test all(round.(fd_M_ssSatCon, digits=3) .== round.(cost_grad_ssSatCon[:, 2:end] , digits=3))

(fd_M_sc_max, fd_M_sc_min, cost_grad_sc) = test_cost_two_spins(GrapeMR.saturation_contrast)
cost_grad_sc_max, cost_grad_sc_min       = cost_grad_sc[1], cost_grad_sc[2]
@test all(round.(fd_M_sc_max, digits=3) .== round.(cost_grad_sc_max[2:end] , digits=3))
@test all(round.(fd_M_sc_min, digits=3) .== round.(cost_grad_sc_min[2:end] , digits=3))

(fd_M_scMx_max, fd_M_scMx_min, cost_grad_scMx) = test_cost_two_spins(GrapeMR.saturation_contrast_Mx)
cost_grad_scMx_max, cost_grad_scMx_min     = cost_grad_scMx[1], cost_grad_scMx[2]
@test all(round.(fd_M_scMx_max, digits=3) .== round.(cost_grad_scMx_max[2:end] , digits=3))
@test all(round.(fd_M_scMx_min, digits=3) .== round.(cost_grad_scMx_min[2:end] , digits=3))

(fd_M_scMt_max, fd_M_scMt_min, cost_grad_scMt) = test_cost_two_spins(GrapeMR.saturation_contrast_Mtrans)
cost_grad_scMt_max, cost_grad_scMt_min     = cost_grad_scMt[1], cost_grad_scMt[2]
@test all(round.(fd_M_scMt_max, digits=3) .== round.(cost_grad_scMt_max[2:end] , digits=3))
@test all(round.(fd_M_scMt_min, digits=3) .== round.(cost_grad_scMt_min[2:end] , digits=3))

