const wandb_project::String = "GrapeMR"
Base.broadcastable(cf::ControlField) = Ref(cf)

"""
    GrapeOutput{T, M1, Mz, F}

Stores the output of GRAPE optimization.

# Fields
- `isochromats::Vector{Isochromat}`: Vector of isochromats containing the spin dynamics.
- `control_field::ControlField{T, M1, Mz}`: Optimized control field after GRAPE optimization.
- `cost_values::Vector{Float64}`: Sequence of cost function values at each iteration.
- `params::Parameters{F}`: Struct containing GRAPE parameters and optimization parameters.
"""
struct GrapeOutput{T<:Real,M1<:AbstractMatrix{T},Mz<:AbstractMatrix{T},F}
    isochromats::Vector{Isochromat}
    control_field::ControlField{T,M1,Mz}
    cost_values::Vector{Float64}
    params::Parameters{F}
end

"""
    grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})

Executes the GRAPE algorithm to optimize control fields for spin dynamics in an NMR/MRI system.

# Arguments
- `p::Parameters`: Struct containing optimization parameters, including maximum iterations and cost function.
- `cf::ControlField`: Initial control field, typically as a spline function.
- `spins::Vector{<:Spins}`: Vector of spins included in the optimization process.

# Returns
- `grape_output::GrapeOutput`: Struct containing optimization results, including optimized control fields, spin dynamics, and cost function values.
"""
function grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    op, gp = p.opt_params, p.grape_params
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter + 1)
    cost_vals = zeros(eltype(cf.B1x), op.max_iter, 1)[:]
    u1x = Matrix{eltype(cf.B1x)}(undef, 1, size(cf.B1x, 2))
    u1y = Matrix{eltype(cf.B1x)}(undef, 1, size(cf.B1x, 2))
    grape_output = GrapeOutput(Vector{Isochromat}(), deepcopy(cf), cost_vals, p)
    ∇x = zeros(eltype(cf.B1x), 1, gp.N)
    ∇y = zeros(eltype(cf.B1x), 1, gp.N)
    mag, adj = zeros(Float64, 4, gp.N + 1), zeros(Float64, 4, gp.N + 1)

    for (ϵ, i) ∈ zip(lr_scheduler, 1:op.max_iter)
        # ϵ   = max(ϵ, eps)
        fill!(∇x, 0.0)
        fill!(∇y, 0.0)

        for spin ∈ spins
            # Forward Propagation 
            iso = dynamics(cf, spin)
            mag = iso.magnetization.dynamics
            cost_vars = gp.cost_function(iso)
            # Cost Variables
            cost, adj_ini = cost_vars
            grape_output.cost_values[i, 1] += cost
            # Adjoint Propagation
            adj = backward_propagation(adj_ini, cf, iso)

            # Save Isochromats from the last iterations
            if i == op.max_iter
                push!(grape_output.isochromats, iso)
            end

            # Gradient
            if gp.fields_opt["B1x"]
                ∇x .+= gradient(adj, mag, Ix)
            end
            if gp.fields_opt["B1y"]
                ∇y .+= gradient(adj, mag, Iy)
            end
        end

        # Control Field
        (u1x, u1y) = update!(cf, (∇x, ∇y), ϵ)
        cf.B1x .= u1x
        cf.B1y .= u1y
    end

    grape_output.control_field.B1x .= u1x
    grape_output.control_field.B1y .= u1y

    # Print Infos
    final_cost = round(grape_output.cost_values[end], digits=3)
    println("\n Final Cost Function Value = $final_cost \n")
    RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


"""
    dynamics(cf::ControlField, spin::Spins)

Calculates and returns the `Isochromat` object representing the spin dynamics under a given control field.

# Arguments
- `cf::ControlField`: Control field affecting the spin dynamics.
- `spin::Spin`: Spin object for which dynamics are computed.

# Returns
- `iso::Isochromat`: Isochromat containing the calculated magnetization dynamics.
"""
function dynamics(cf::ControlField, spin::Spins)
    mag = forward_propagation(cf, spin)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spin)
    return iso
end


"""
    gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)

Calculates the gradient of the cost function with respect to the Hamiltonian for each time step.

# Arguments
- `χ::Matrix{Float64}`: Adjoint state matrix.
- `M::Matrix{Float64}`: Forward propagation matrix.
- `H::Matrix`: Hamiltonian matrix.

# Returns
- `grad::Matrix{Float64}`: Gradient of the cost function, as a 1xN matrix.
"""
function gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::AbstractMatrix{Int64})
    # TODO: refactor this as gradient!(grad, χ, M, H)
    grad = zeros(Float64, 1, size(M, 2) - 1)
    for i ∈ 1:(size(M, 2)-1)
        grad[1, i] = dot(
            transpose(view(χ, :, i + 1)),
            H,
            view(M, :, i + 1)
        )
    end
    return grad
end

"""
    update!(cf::ControlField, ∇xy::Tuple, ϵ::Float64)

Updates the control fields based on the calculated gradient and a learning rate.

# Arguments
- `cf::ControlField`: Control field struct to be updated.
- `∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}`: Gradients for the x and y components of the field.
- `ϵ::Float64`: Learning rate for gradient descent.

# Returns
- `(u1x, u1y)`: Updated x and y control fields.
"""
function update!(cf::ControlField, ∇xy::Tuple{Matrix{Float64},Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ .* ∇xy[1]
    u1y = cf.B1y .- ϵ .* ∇xy[2]
    return u1x, u1y
end

"""
    run_grape_optimization(config_path::String)

Runs the GRAPE optimization process based on a TOML configuration file. Initializes spins, sets up parameters, executes optimization, and optionally saves and plots the results.

# Arguments
- `config_path::String`: Path to the TOML configuration file containing parameters for spins, optimization, and control fields.

# Configuration File Structure
The TOML configuration file should include sections like:
    - **spins**: Defines spin properties.
    - **grape_parameters**: Parameters for the GRAPE optimization.
    - **optimization_parameters**: Parameters for the hyperparameter optimization.
    - **control_field**: Control field specifications.
    - **save_files**: Settings for saving outputs.
    - **plot**: Plot settings.

# Returns
- Produces and saves results depending on configuration settings, including optimized control fields, cost values, optional Bruker export, and plots.

# Example
```julia
run_grape_optimization("path/to/config.toml")
"""
function run_grape_optimization(config_path::String)
    tm = TOML.parsefile(config_path)
    @info "Configuration:"
    pprintln(tm)

    offsets = collect(-tm["spins"]["offset"]:1:tm["spins"]["offset"])

    # Spin Object
    spins = GrapeMR.Spin(
        tm["spins"]["M0"],
        [s["T1"] for s in tm["spins"]["intrinsics"]],
        [s["T2"] for s in tm["spins"]["intrinsics"]],
        offsets, tm["spins"]["delta_B1"],
        [s["target"] for s in tm["spins"]["intrinsics"]],
        [s["label"] for s in tm["spins"]["intrinsics"]]
    )

    # Grape Parameters
    mask_dict = Dict(k => Bool(v) for (k, v) ∈ tm["grape_parameters"]["fields2optimize"])
    grape_params = GrapeParams(
        tm["grape_parameters"]["time_steps"],
        eval(Symbol(tm["grape_parameters"]["cost_function"])),
        mask_dict
    )

    # Optimization Parameters
    if tm["optimization_parameters"]["hyper_opt"]
        hyper_opt = bohb_hyperopt(spins, grape_params, LinRange(0.01, 0.5, 9), 2187)
        # hyper_opt = random_hyperopt(spins, grape_params, LinRange(0.01, 1.0, 15), range(500, 2000, step = 100))
        Tc, poly_start, poly_degree, max_iter = hyper_opt.minimizer
        opt_params = OptimizationParams(
            poly_start,
            poly_degree,
            Int(ceil(max_iter))
        )
        # Initial RF Pulse Object
        control_field = spline_RF(grape_params.N, Tc, tm["control_field"]["B1ref"])
    else
        opt_params = OptimizationParams(
            tm["optimization_parameters"]["poly_start"],
            tm["optimization_parameters"]["poly_degree"],
            Int(ceil(tm["optimization_parameters"]["max_iter"]))
        )
        # Initial RF Pulse Object
        control_field = spline_RF(grape_params.N, tm["control_field"]["control_time"], tm["control_field"]["B1ref"])
    end

    # Parameters 
    params = Parameters(grape_params, opt_params)

    # Run Optimization
    grape_output = grape(params, control_field, spins)

    # Save data
    if tm["save_files"]["enabled"]
        # Save output data
        if tm["optimization_parameters"]["hyper_opt"]
            experiment_folder = save_grape_data(grape_output; folder_path=tm["save_files"]["folder_path"])
            experiment_folder = save_hyperopt_data(hyper_opt; folder_path=tm["save_files"]["folder_path"])
        else
            experiment_folder = save_grape_data(grape_output; folder_path=tm["save_files"]["folder_path"])
        end
        # Export Bruker data
        if tm["save_files"]["export_bruker"]
            export_bruker(grape_output; folder_path=tm["save_files"]["bruker_folder_path"])
        end
    end

    if tm["plot"]
        # Plots
        display(plot_cost_values(grape_output.cost_values, grape_params))
        display(plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats))
        display(plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control))
        # TODO add if s[t1] > 2 plot(iso[end]) to get time dynamics of the second spin
    end
end
