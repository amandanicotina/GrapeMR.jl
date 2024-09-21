using Test
using GrapeMR
using FiniteDifferences
using Plots

# function test_gradient_Bx()
    # Parameters
    N   = 1000;
    t_c = 0.5; #[s]
    cost_func = :target_one_spin
    grape_params = GrapeParams(N, cost_func, [true true false])

    # RF
    B1ref = 1.0
    B1x   = spline_RF(grape_params.N, t_c)'
    B1y   = spline_RF(grape_params.N, t_c)'
    Bz    = zeros(1, grape_params.N)

    # Spin
    M0  = [0.0, 0.0, 1.0];
    T1  = [0.8, 0.4];
    T2  = [0.6, 0.3];
    B0  = [0.0];
    ΔB1 = [1.0];
    target = ["min", "max"];
    label  = ["T1-$(Int(T1[1]*1e3))ms", "T1-$(Int(T1[2]*1e3))ms"]

    # Spin and RF objects
    cf_tmp = ControlField(B1x, B1y, B1ref, Bz, t_c)
    spins  = GrapeMR.Spin(M0, T1, T2, B0, ΔB1, target, label)

    mag_vals = forward_propagation(cf_tmp, spins[1])
    dyn_vals = Magnetization(mag_vals)
    iso_vals = Isochromat(dyn_vals, spins[1])
    cost_val, cost_grad = GrapeMR.cost_function(iso_vals, grape_params.cost_function) 
    adj = backward_propagation(cf_tmp, iso_vals, cost_grad)

    # True Gradient
    true_grad = GrapeMR.gradient(adj, mag_vals, Iy)

    # Finite difference
    Δcf = 1e-2
    fd_cf = finite_difference_field_symmetric(cost_func, cf_tmp, spins[1], Δcf)

    # f(x) = cost_function(x, cost_func)
    # grad(central_fdm(5, 1), f, iso_vals)

    plot(fd_cf')
    plot!(true_grad')
    
# end

# using FiniteDifferences
# central_diff = FiniteDifferences.central_fdm(2, 1)
# cf =  B1x -> cost_wrapper(B1x)[1]
# fd_cf = FiniteDifferences.grad(central_diff, cf, B1x)  







