ctx = GrapeContext(
    control_field,
    spins,
    fields_opt,
    grape_params,
    zeros(T, 4, n+1),   # magnetization
    zeros(T, 4, n+1)    # adjoint
)


get_max_iter(config::GradientDescentConfig) = config.max_iter
get_max_iter(config::ManualGradientDescentConfig) = config.max_iter
get_max_iter(config::BFGSConfig) = config.max_iter


function cost_function_vec!(u_vec, control_field, spins, magnetization, cost_function, fields_opt)
    n = size(control_field.B1x, 2)
    if fields_opt["B1x"]
        control_field.B1x .= reshape(view(u_vec, 1:n), 1, :)
    end
    if fields_opt["B1y"]
        control_field.B1y .= reshape(view(u_vec, n+1:2n), 1, :)
    end

    total_cost = 0.0
    for spin in spins
        forward_propagation!(magnetization, control_field, spin)
        iso = Isochromat(Magnetization(magnetization), spin)
        cost, _ = cost_function(iso)
        total_cost += cost
    end
    return total_cost
end

function gradient_function_vec!(G, u_vec, control_field, spins, magnetization, adjoint, cost_function, fields_opt)
    n = size(control_field.B1x, 2)
    if fields_opt["B1x"]
        control_field.B1x .= reshape(view(u_vec, 1:n), 1, :)
    end
    if fields_opt["B1y"]
        control_field.B1y .= reshape(view(u_vec, n+1:2n), 1, :)
    end

    grad_x = zeros(Float64, 1, n)
    grad_y = zeros(Float64, 1, n)

    for spin in spins
        forward_propagation!(magnetization, control_field, spin)
        iso = Isochromat(Magnetization(magnetization), spin)
        _, adj_init = cost_function(iso)
        backward_propagation!(adjoint, control_field, iso, adj_init)

        if fields_opt["B1x"]
            grad_x .+= gradient(adjoint, magnetization, Ix)
        end
        if fields_opt["B1y"]
            grad_y .+= gradient(adjoint, magnetization, Iy)
        end
    end

    G .= vcat(vec(grad_x), vec(grad_y))
    return G
end