using Test, SafeTestsets
using GrapeMR

function test_cost_one_spin(cost_func::Symbol; target::String = "-")
    # Parameters
    N = 1000;
    grape_params = GrapeParams(N, cost_func, [true true false])

    # RFs
    t_c   = 0.5; #[s]
    B1ref = 1.0

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8];
    T2  = [0.6];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = [target];
    label  = ["T1-$(Int(T1[1]*1e3))ms"]

    # Spin and RF objects
    control_field = spline_RF(grape_params.N, t_c, B1ref)
    spins   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    iso_vec = dynamics.(control_field, spins)

    # Cost gradient
    cost_vars = GrapeMR.cost_function.(iso_vec, grape_params.cost_function)
    cost_grad = getindex.(cost_vars, 2)
   
    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(grape_params.cost_function, iso_vec[1], ΔM)

    return fd_M, cost_grad
end

function test_cost_two_spins(cost_func::Symbol)
    # Parameters
    N = 1000;
    grape_params = GrapeParams(N, cost_func, [true true false])

    # RFs
    t_c   = 0.5; #[s]
    B1ref = 1.0

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8, 0.4];
    T2  = [0.6, 0.3];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = ["max", "min"];
    label  = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

    # Spin and RF objects
    control_field = spline_RF(grape_params.N, t_c, B1ref)
    spins   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    iso = dynamics.(Ref(control_field), spins)

    # Cost gradient
    cost_vars = GrapeMR.cost_function.(iso, grape_params.cost_function)
    cost_grad = getindex.(cost_vars, 2)

    # Finite difference
    ΔM   = 1e-10;
    fd_M_max = finite_difference_cost(grape_params.cost_function, iso[1], ΔM)
    fd_M_min = finite_difference_cost(grape_params.cost_function, iso[2], ΔM)

    return fd_M_max, fd_M_min, cost_grad
end

(fd_M_en, cost_grad_en) = test_cost_one_spin(:euclidean_norm)
@test all(round.(fd_M_en, digits=3) .== round.(cost_grad_en[:, 2:end], digits=3))

(fd_M_tos, cost_grad_tos) = test_cost_one_spin(:spin_target; target = "[0.0, 1.0, 0.0]")
@test all(round.(fd_M_tos, digits=3) .== round.(cost_grad_tos[:, 2:end] , digits=3))

(fd_M_sc_max, fd_M_sc_min, cost_grad_sc) = test_cost_two_spins(:saturation_contrast)
cost_grad_sc_max, cost_grad_sc_min       = cost_grad_sc[1], cost_grad_sc[2]
@test all(round.(fd_M_sc_max, digits=3) .== round.(cost_grad_sc_max[2:end] , digits=3))
@test all(round.(fd_M_sc_min, digits=3) .== round.(cost_grad_sc_min[2:end] , digits=3))

(fd_M_scMx_max, fd_M_scMx_min, cost_grad_scMx) = test_cost_two_spins(:saturation_contrast_Mx)
cost_grad_scMx_max, cost_grad_scMx_min     = cost_grad_scMx[1], cost_grad_scMx[2]
@test all(round.(fd_M_scMx_max, digits=3) .== round.(cost_grad_scMx_max[2:end] , digits=3))
@test all(round.(fd_M_scMx_min, digits=3) .== round.(cost_grad_scMx_min[2:end] , digits=3))

(fd_M_scMt_max, fd_M_scMt_min, cost_grad_scMt) = test_cost_two_spins(:saturation_contrast_Mtrans)
cost_grad_scMt_max, cost_grad_scMt_min     = cost_grad_scMt[1], cost_grad_scMt[2]
@test all(round.(fd_M_scMt_max, digits=3) .== round.(cost_grad_scMt_max[2:end] , digits=3))
@test all(round.(fd_M_scMt_min, digits=3) .== round.(cost_grad_scMt_min[2:end] , digits=3))