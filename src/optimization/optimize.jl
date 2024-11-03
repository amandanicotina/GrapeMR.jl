const wandb_project::String = "GrapeMR"
Base.broadcastable(cf::ControlField) = Ref(cf)


struct GrapeOutput
    isochromats::Vector{Isochromat}
    control_field::ControlField
    cost_values::Union{Array{Float64}, Metal.MtlVector{Float32, Metal.SharedStorage}}
    params::Parameters
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
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter+1) 
    cost_vals    = zeros(Float64, op.max_iter, 1)[:]
    u1x, u1y     = zeros(Float64, length(cf.B1x)), zeros(Float64, length(cf.B1x))
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
    ∇x  = zeros(Float64, 1, gp.N)
    ∇y  = zeros(Float64, 1, gp.N)
    mag, adj = zeros(Float64, 4, gp.N+1), zeros(Float64, 4, gp.N+1)

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
            adj = backward_propagation(cf, iso, adj_ini)
            
            # Save Isochromats from the last iterations
            if i == op.max_iter
                push!(grape_output.isochromats, iso)
            end
            
            # Gradient
            if gp.fields_opt[1]
                ∇x .+= gradient(adj, mag, Ix)
            end 
            if gp.fields_opt[2]
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
    final_cost = round(grape_output.cost_values[end], digits = 3)
    # println("\n Final Cost Function Value = $final_cost \n")
    # RF_pulse_analysis(grape_output.control_field)

    return grape_output
end


function metal_grape(p::Parameters, cf::ControlField, spins::Vector{<:Spins})
    # Extract parameters
    op, gp = p.opt_params, p.grape_params
    lr_scheduler = Poly(start=op.poly_start, degree=op.poly_degree, max_iter=op.max_iter+1)
    
    # Initialize MetalBuffers for unified memory
    cost_vals = Metal.zeros(Float32, op.max_iter; storage = Metal.SharedStorage)  # Cost values shared between CPU and GPU
    u1x = Metal.zeros(Float32, length(cf.B1x); storage = Metal.SharedStorage)
    u1y = Metal.zeros(Float32, length(cf.B1x); storage = Metal.SharedStorage)
    grape_output = GrapeOutput([], deepcopy(cf), cost_vals, p)
    ∇x = Metal.zeros(Float32, 1, gp.N; storage = Metal.SharedStorage)
    ∇y = Metal.zeros(Float32, 1, gp.N; storage = Metal.SharedStorage)
    adj = Metal.zeros(Float32, 4, gp.N+1; storage = Metal.SharedStorage)

    # Main optimization loop
    for (ϵ, i) ∈ zip(lr_scheduler, 1:op.max_iter)
        Metal.fill!(∇x, 0.0f0)
        Metal.fill!(∇y, 0.0f0)

        isos = Metal.zeros(Float32, length(spins), 4, gp.N+1; storage = Metal.SharedStorage)
        # TODO: Ideally we would have a function that does this in a single kernel
        for (idx, spin) ∈ enumerate(spins)
            iso = isos[idx, :, :]
            # Forward Propagation (ensure `dynamics` and `cost_function` are compatible with Metal)
            metal_dynamics!(cf, spin, iso)  # may need Metal compatibility
            # Assumption that needs to be tested: this is fast becuase of the unified memory
            cpu_iso = unsafe_wrap(Matrix{Float32}, iso, size(iso))
            cost, adj_ini = gp.cost_function(Isochromat(Magnetization(Float64.(cpu_iso)), spin))
            adj_ini = Float32.(adj_ini)
            grape_output.cost_values[i] += cost

            # Adjoint Propagation (compatible with Metal)
            χ = MtlMatrix{Float32, Metal.SharedStorage}(undef, 4, gp.N+1)
            χ .= adj_ini
            metal_backward_propagation(
                Float32.(cf.t_control), Float32.(cf.B1x), Float32.(cf.B1y), Float32.(cf.Bz), Float32.(spin.B0inho), Float32.(spin.B1inho), Float32.(spin.T1), Float32.(spin.T2) , χ
            )

            # if i == op.max_iter
            #     push!(grape_output.isochromats, iso)
            # end

            # Gradient calculation in Metal (ensure `gradient` is Metal-compatible)
            # grad = Metal.zeros(Float32, 1, length(iso[1,:])-1; storage=Metal.SharedStorage)
            # if gp.fields_opt[1]
            #     ∇x .+= metal_gradient(χ, iso, metal_Ix, grad)
            # end
            # if gp.fields_opt[2]
            #     ∇y .+= metal_gradient(χ, iso, metal_Iy, grad)
            # end
        end

        # Update control fields
        (u1x, u1y) = metal_update!(Float32.(cf.B1x), Float32.(cf.B1y), (∇x, ∇y), Float32.(ϵ))
        cf.B1x .= u1x
        cf.B1y .= u1y
    end

    # Finalize control field values
    grape_output.control_field.B1x .= Array(u1x)  # Convert back to standard arrays if necessary
    grape_output.control_field.B1y .= Array(u1y)

    # Print final cost for reference
    final_cost = round(grape_output.cost_values[end], digits = 3)

    return grape_output
end


"""
    metal_dynamics!(cf::ControlField, spins::Spin)

Function that returns the Isochromat object with the already calculated dynamics.

# Input
- cf::ControlField - Adjoint State
- spins::Spin

# Output
- nothing

Note: Metal.jl kernels are not allowed to return values, so we need to pass the magnetization matrix as an argument
"""
function metal_dynamics!(cf::ControlField, spin::Spin, M::MtlMatrix{Float32, Metal.SharedStorage})
    metal_forward_propagation(
        Float32.(cf.t_control), Float32.(cf.B1x), Float32.(cf.B1y), Float32.(cf.Bz), Float32.(spin.B0inho), Float32.(spin.B1inho), Float32.(spin.T1), Float32.(spin.T2), Float32.(spin.M_init),
        M
    )
    return
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
function dynamics(cf::ControlField, spin::Spin)
    mag = forward_propagation(cf, spin)
    dyn = GrapeMR.Magnetization(mag)
    iso = Isochromat(dyn, spin)
    return iso
end


"""
    gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)

# Arguments  
- χ = (::Matrix{Float64}) - Adjoint State
- M = (::Matrix{Float64}) - Forward Propagation
- H = (::Matrix) - Hamiltonian

# Outputs
- ΔJ - 1xN matrix
"""
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

function metal_gradient(
    χ::MtlMatrix{Float32, Metal.SharedStorage}, 
    M::MtlMatrix{Float32, Metal.SharedStorage}, 
    H::MtlMatrix{Float32, Metal.SharedStorage},
    grad::MtlMatrix{Float32,Metal.SharedStorage}
    )
    for i ∈ 1:(length(M[1,:])-1)
        grad[1,i] = transpose(χ[:,i+1])*H*M[:,i+1]
    end
    return
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
function update!(cf::ControlField, ∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ*∇xy[1]
    u1y = cf.B1y .- ϵ*∇xy[2]
    return u1x, u1y
end

function metal_update!(cf_B1x::Matrix{Float32}, cf_B1y::Matrix{Float32}, ∇xy::Tuple{MtlArray{Float32, 2, Metal.SharedStorage}, MtlArray{Float32, 2, Metal.SharedStorage}}, ϵ::Float32)
    u1x = cf_B1x .- ϵ*∇xy[1]
    u1y = cf_B1y .- ϵ*∇xy[2]
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

    offsets = collect(-tm["spins"]["offset"]:5:tm["spins"]["offset"])

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
    mask = tm["grape_parameters"]["fields2optimize"]
    grape_params = GrapeParams(
                    tm["grape_parameters"]["time_steps"], 
                    Symbol(tm["grape_parameters"]["cost_function"]), 
                    reshape(mask, 1, length(mask))  # we need a 1x3 Matrix{Bool} instead of a Vector{Bool}
                    )

    # Optimization Parameters
    if tm["optimization_parameters"]["hyper_opt"]
        hyper_opt = bohb_hyperopt(spins, grape_params, LinRange(0.01, 0.1, 9), 2187)
        # hyper_opt = random_hyperopt(spins, grape_params, LinRange(0.01, 1.0, 15), range(500, 2000, step = 100))
        Tc, poly_start, poly_degree, max_iter = hyper_opt.minimizer 
        opt_params   = OptimizationParams(
                        poly_start, 
                        poly_degree, 
                        Int(ceil(max_iter))
                        )
        # Initial RF Pulse Object
        control_field = spline_RF(grape_params.N, Tc, tm["control_field"]["B1ref"]) 
    else
        opt_params   = OptimizationParams(
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
    grape_output = @time grape(params, control_field, spins); 

    # Save data
    if tm["save_files"]["enabled"]
        # Save output data
        if tm["optimization_parameters"]["hyper_opt"]
            experiment_folder = save_grape_data(grape_output; folder_path = tm["save_files"]["folder_path"])
            experiment_folder = save_hyperopt_data(hyper_opt; folder_path = tm["save_files"]["folder_path"])
        else
            experiment_folder = save_grape_data(grape_output; folder_path = tm["save_files"]["folder_path"])
        end
        # Export Bruker data
        if tm["save_files"]["export_bruker"]
            export_bruker(grape_output; folder_path = tm["save_files"]["bruker_folder_path"])
        end
    end

    if tm["plot"]
        # Plots
        plot_cost_values(grape_output.cost_values, grape_params)
        plot_magnetization_2D(grape_output.isochromats)
        plot_magnetization_control_field(grape_output.control_field, grape_output.isochromats)
        plot_control_fields(grape_output.control_field; unit="Hz")
        plot_magnetization_time(grape_output.isochromats[1], grape_output.control_field.t_control)
        plot_transverse_time(grape_output.isochromats, grape_output.control_field.t_control)
        plot_transverse_magnetization(grape_output.isochromats)
    end


    # go = load_grape_data(experiment_folder)
    spins = GrapeMR.Spin(
                    tm["spins"]["M0"],
                    [1e8, 1e8], 
                    [1e8, 1e8],
                    offsets, tm["spins"]["delta_B1"], 
                    [s["target"] for s in tm["spins"]["intrinsics"]],
                    [s["label"] for s in tm["spins"]["intrinsics"]]
                )

    cf = grape_output.control_field

    iso_noRelax = dynamics.(cf, spins)    
    plot_magnetization_control_field(grape_output.control_field, iso_noRelax)
    plot_transverse_time(iso_noRelax, grape_output.control_field.t_control)
    plot_transverse_magnetization(iso_noRelax)

end
