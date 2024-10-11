using Test
using GrapeMR

function test_cost_one_spin(cost_func::Symbol)
    # Parameters
    N   = 1000;
    t_c = 0.5; #[s]
    grape_params = GrapeParams(N, cost_func, [true true false])

    # RFs
    B1ref = 1.0
    B1x = spline_RF(grape_params.N, t_c)'
    B1y = spline_RF(grape_params.N, t_c)'
    Bz  = zeros(1, grape_params.N)
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8];
    T2  = [0.6];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = ["-"];
    label  = ["T1-$(Int(T1[1]*1e3))ms"]

    # Spin and RF objects
    control_field = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag = forward_propagation(control_field, spins[1])
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spins[1])

    # Cost gradient
    (_, cost_grad) = GrapeMR.cost_function(iso, grape_params.cost_function)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(grape_params.cost_function, iso, ΔM)

    return fd_M, cost_grad
end


function test_cost_two_spins(cost_func::Symbol)
    # Parameters
    N   = 1000;
    t_c = 0.5; #[s]
    grape_params = GrapeParams(N, cost_func, [true true false])

    # RFs
    B1ref = 1.0
    B1x = spline_RF(grape_params.N, t_c)'
    B1y = spline_RF(grape_params.N, t_c)'
    Bz  = zeros(1, grape_params.N)
    control_field = ControlField(B1x, B1y, B1ref, Bz, t_c)

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8, 0.4];
    T2  = [0.6, 0.3];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = ["min", "max"];
    label  = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

    # Spin and RF objects
    control_field = ControlField(B1x, B1y, 1.0, Bz, t_c)
    spins   = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag = forward_propagation(control_field, spins[1])
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spins[1])

    # Cost gradient
    (_, cost_grad) = GrapeMR.cost_function(iso, grape_params.cost_function)

    # Finite difference
    ΔM   = 1e-10;
    fd_M = finite_difference_cost(grape_params.cost_function, iso, ΔM)

    return fd_M, cost_grad
end


(fd_M_en, cost_grad_en) = test_cost_one_spin(:euclidean_norm)
@test all(round.(fd_M_en, digits=3) .== round.(cost_grad_en[2:end,:], digits=3))

(fd_M_tos, cost_grad_tos) = test_cost_one_spin(:spin_target)
@test all(round.(fd_M_tos, digits=3) .== round.(cost_grad_tos[2:end,:] , digits=3))

(fd_M_sc, cost_grad_sc) = test_cost_two_spins(:saturation_contrast)
@test all(round.(fd_M_sc, digits=3) .== round.(cost_grad_sc[2:end,:] , digits=3))

(fd_M_scMx, cost_grad_scMx) = test_cost_two_spins(:saturation_contrast_Mx)
@test all(round.(fd_M_scMx, digits=3) .== round.(cost_grad_scMx[2:end,:] , digits=3))

(fd_M_scMt, cost_grad_scMt) = test_cost_two_spins(:saturation_contrast_Mtrans)
@test all(round.(fd_M_scMt, digits=3) .== round.(cost_grad_scMt[2:end,:] , digits=3))