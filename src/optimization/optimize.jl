const wandb_project::String = "GrapeMR"
Base.broadcastable(cf::ControlField) = Ref(cf)

struct GrapeOutput{T<:Real, M1<:AbstractMatrix{T}, Mz<:AbstractMatrix{T}, F}
    isochromats::Vector{Isochromat}
    control_field::ControlField{T, M1, Mz}
    cost_values::Vector{Float64}
    params::Parameters{F}
end



"""
    grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})

Implements the grape algorithm [khaneja2005optimal](@cite) 

# Arguments
- `op::OptimizationParams`: Parameters for the optimization itself: max iterations, 
- `gp::GrapeParams`: Parameters related to Grape itself: time points, cost function, mask for which fields are being optimized.
- `cf::ControlField`: Initial control field - spline function -
- `spins::Vector{<:Spins}`: Vector with all spins included in the optimization

# Outputs
A scruct cointaing all optimization results:
- `grape_output::GrapeOutput': Data type with the optimized control fields, spin information and spin dynamics.
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
            grape_output.cost_values[i,1] += cost
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
    # println("\n Final Cost Function Value = $final_cost \n")
    # RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


"""
    dynamics(cf::ControlField, spins::Spin)

Function that returns the Isochromat object with the already calculated dynamics.

    # Input  
    - cf::ControlField - Adjoint State
    - spins::Spin

    # Output
    - iso::Isochromat
"""
function dynamics(cf::ControlField, spin::Spins)
    mag = forward_propagation(cf, spin)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spin)
    return iso
end


"""
    gradient!(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)

# Arguments  
- χ = (::Matrix{Float64}) - Adjoint State
- M = (::Matrix{Float64}) - Forward Propagation
- H = (::Matrix) - Hamiltonian

# Outputs
- ΔJ - 1xN matrix
"""
# function gradient!(grad::AbstractMatrix, χ::AbstractMatrix, M::AbstractMatrix, H::AbstractMatrix{Int64})
#     for i ∈ 1:(size(M, 2)-1)
#         grad[1, i] = dot(
#             transpose(view(χ, :, i + 1)),
#             H,
#             view(M, :, i + 1)
#         )
#     end
#     return grad
# end
function gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::AbstractMatrix{Int64})
    # TODO: refactor this as gradient!(grad, χ, M, H)
    grad = zeros(Float64, 1, size(M, 2)-1)
    for i ∈ 1:(size(M, 2)-1)
        grad[1,i] = dot(
            transpose(view(χ, :,i+1)), 
            H,
            view(M, :,i+1)
        )
    end
    return grad
end

"""
    update(cf::ControlField, ∇xy::Tuple, ϵ::Float64)

update
# Arguments  
- cf:  (::ControlField) - Control fields struct
- ∇xy: (::Tuple) - Calculated gradients for x and y components
- ϵ:   (::Float64) - Weigth of gradient

# Outputs
- Control Field - 1xN matrix
"""
function update!(cf::ControlField, ∇xy::Tuple{Matrix{Float64},Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ .* ∇xy[1]
    u1y = cf.B1y .- ϵ .* ∇xy[2]
    return u1x, u1y
end

"""
    run_grape_optimization(config_path::String)

This function runs the GRAPE optimization process for an NMR/MRI system based on a given configuration file. The function initializes spin systems, sets up optimization parameters, and executes the optimization using the GRAPE algorithm. The results are optionally saved and plotted based on the configuration file.

# Arguments
- `config_path::String`: 
    Path to the TOML configuration file. This file should contain all necessary parameters for spins, optimization, and control fields.

# Configuration File Structure
    The configuration file should follow the TOML format and include the following sections:

    - **spins**:
    - **grape_parameters**:
    - **optimization_parameters**:
    - **control_field**:
    - **save_files**:
    - **plot**:

# Outputs
    The function produces several outputs depending on the configuration:
    - GRAPE optimization results, including optimized control fields and cost values.
    - Hyperparameter optimization results, if enabled.
    - Optional saving of output data and Bruker export.
    - Plots visualizing the magnetization and control fields over time.

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
    mask_dict = Dict(k => Bool(v) for (k,v) ∈ tm["grape_parameters"]["fields2optimize"])
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
    grape_output = @time grape_const_epsilon(params, control_field, spins)

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
        plot_magnetization_2D(grape_output.isochromats)
        display(plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats))
        plot_control_fields(grape_output.control_field; unit="Hz")
        plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control)
        # plot_transverse_time(grape_output.isochromats, grape_output.control_field.t_control)
        # plot_transverse_magnetization(grape_output.isochromats)
    end


    # spins = GrapeMR.Spin(
    #     tm["spins"]["M0"],
    #     [1e8, 1e8],
    #     [1e8, 1e8],
    #     offsets, tm["spins"]["delta_B1"],
    #     [s["target"] for s in tm["spins"]["intrinsics"]],
    #     [s["label"] for s in tm["spins"]["intrinsics"]]
    # )

    # cf = grape_output.control_field

    # iso_noRelax = dynamics.(cf, spins)
    # plot_magnetization_control_field(grape_output.control_field, iso_noRelax)
    # plot_transverse_time(iso_noRelax, grape_output.control_field.t_control)
    # plot_transverse_magnetization(iso_noRelax)

end
