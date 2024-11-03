"""
    calculate_cost_gradient(cost_func::Num, dict::Dict, vars::Vector{Num})

Generic function that maps variables into the gradient calculation using the Symbolics package
    dM = Differential(vars[])
    dCost_M = expand_derivatives(dM(cost_func))
    P = Symbolics.unwrap(substitute.(dCost_M, (dict,))[1])

# Arguments
- 'cost_func::Num': Cost Function from pre-defined dictionary
- 'dict::Dict': Variable dictionary
- 'vars::Vector{Num}': Function variables, Mx, My, etc.

# Outputs
- 'cost_gradient::Vector{Float64}': 
"""
function calculate_cost_gradient(cost_func::Num, dict::Dict, vars::Vector{Num})
    error("deleted")
end


"""
    get_cost_and_gradient(iso::Isochromat, cost_expr::Num, vars::Vector{Num})

Generic function that maps variables into the gradient calculation using the Symbolics package
    dM = Differential(vars[1])
    dCost_M = expand_derivatives(dM(cost_func))
    P = Symbolics.unwrap(substitute.(dCost_M, (dict,))[1])

# Arguments
- 'iso::Isochromat': Isochromts
- 'cost_expr::Num': Cost Function
- 'vars::Vector{Num}': Function variables -> Mx, My, etc.

# Outputs
- 'cost_func_val::Float64': Cost Function Value
- 'cost_gradient::Vector{Float64}': 
"""
function get_cost_and_gradient(iso::Isochromat, cost_expr::Num, vars::Vector{Num})
    # TODO: replace with ForwardDiff on a StaticArrays.SVector
    error("to replace with ForwardDiff")
end 

############################################################################################
#                                     Cost Functions                                       #
############################################################################################

euclidean_norm_cost(x::AbstractVector, Nspins) = sqrt(sum(abs2, x) + 1e-15) / Nspins

function euclidean_norm(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    x = SVector(m[2, end], m[3, end], m[4, end])
    val = euclidean_norm_cost(x, s.Nspins)
    grad = ForwardDiff.gradient(Base.Fix2(euclidean_norm_cost, s.Nspins), x)
    return val, vcat(zero(eltype(grad)), grad)
end

function cost_func_do_dawid(iso::Isochromat)
    vars = @variables Mx, My, Mz
    cost = (sin(Mx) - cos(My))*atan(Mz)
    return get_cost_and_gradients(iso, cost, vars)
end

function spin_target(iso::Isochromat; M_tar = [0.0, 1.0, 0.0])
    s = iso.spin
    vars  = @variables Mx, My, Mz
    M_tar = eval(Meta.parse(iso.spin.target))  # TODO: remove!!!! danger!!!!
    Mx_tar, My_tar, Mz_tar = M_tar
    cost_expr = sqrt((Mx - Mx_tar)^2 + (My - My_tar)^2 + (Mz - Mz_tar)^2 + 1e-15) / s.Nspins
    return get_cost_and_gradient(iso, cost_expr, vars)
end

function saturation_contrast(iso::Isochromat)
    vars = @variables Mx, My, Mz
    s = iso.spin
    if s.target == "max"
        cost_expr = (1 - sqrt(Mz^2 + 1e-15)) / s.Nspins
    elseif s.target == "min"
        cost_expr = sqrt(Mx^2 + My^2 + Mz^2 + 1e-15) / s.Nspins
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return get_cost_and_gradient(iso, cost_expr, vars)
end

function saturation_contrast_Mx(iso::Isochromat)
    vars = @variables Mx, My, Mz
    s = iso.spin
    if s.target == "max"
        cost_expr = (1 - sqrt(Mx^2 + 1e-15)) / s.Nspins 
    elseif s.target == "min"
        cost_expr = sqrt(Mx^2 + My^2 + 1e-15) / s.Nspins
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return get_cost_and_gradient(iso, cost_expr, vars)
end

function saturation_contrast_Mtrans(iso::Isochromat)
    vars = @variables Mx, My, Mz
    s = iso.spin

    if s.target == "max"
        cost_expr = (1 - sqrt(Mx^2 + My^2 + 1e-15)) / s.Nspins
    elseif s.target == "min"
        cost_expr = sqrt(Mx^2 + My^2 + 1e-15) / s.Nspins
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return get_cost_and_gradient(iso, cost_expr, vars)
end

function steady_state_offset_targets(iso::Isochromat) 
    vars = @variables Mx, My, Mz
    s = iso.spin

    ss = steady_state_matrix(iso)
    Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

    cost_expr = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)/s.Nspins

    return get_cost_and_gradient(iso, cost_expr, vars)
end

function saturation_contrast_steady_state(iso::Isochromat)
    vars = @variables Mx, My, Mz
    s = iso.spin
    if s.target == "min"
        cost_expr = sqrt(Mx^2 + My^2 + Mz^2 + 1e-15) / s.Nspins
    elseif s.target == "max"
        # Steady State
        ss = steady_state_matrix(iso)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

        # cost_expr = (1 - sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15))/s.Nspins
        cost_expr = (1 - sqrt((Mz - Mz_ss)^2 + 1e-15))/s.Nspins
    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end

    return get_cost_and_gradient(iso, cost_expr, vars)
end

function cos_func_do_dawid(iso::Isochromat)
    vars = @variables Mx, My, Mz
    cost = (sin(Mx) - cos(My))*atan(Mz)
    return get_cost_and_gradients(iso, cost, vars)
end


"""
Cost Function's dictionary
"""
const COST_FUNCTIONS = Dict(
    :euclidean_norm => euclidean_norm,
    :spin_target => spin_target,
    :steady_state_offset_targets => steady_state_offset_targets,
    :saturation_contrast => saturation_contrast,
    :saturation_contrast_Mx => saturation_contrast_Mx,
    :saturation_contrast_Mtrans => saturation_contrast_Mtrans,
    :saturation_contrast_steady_state => saturation_contrast_steady_state
)




# function saturation_contrast_square(iso::Isochromat)
#     c = 0.0;
#     m = iso.magnetization.dynamics
#     s = iso.spin

#     if s.target == "max"
#         c = (1 - sum(m[4,end]*m[4,end]))/s.Nspins
#     elseif s.target == "min"
#         c = sum(m[2:end,end].*m[2:end,end])/s.Nspins
#     end
#     return c
# end

# function target_steady_state(iso::Isochromat) 
#     c = 0
#     s = iso.spin
#     if s.target == "min"
#         # Steady State
#         ss = steady_state_matrix(s)
#         Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
#         # Magnetization
#         mag = iso.magnetization.dynamics
#         Mx  = mag[2,end]
#         My  = mag[3,end]
#         Mz  = mag[4,end]
#         c = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)

#     elseif s.target == "max"
#         # Steady State
#         ss = steady_state_matrix(s)
#         Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
#         # Magnetization
#         mag = iso.magnetization.dynamics
#         Mx  = mag[2,end]
#         My  = mag[3,end]
#         Mz  = mag[4,end]
#         c = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)

#     else
#         error(" $(s.target) is not a matching target. Valid targets are max or min")
#     end

#     return c/s.Nspins
# end



# # """
# # Old Functions
# # """
# # function target_phase_encoding(iso::Isochromat; M_tar = [0.0, 1.0, 0.0])
# #     s  = iso.spin
# #     m  = iso.magnetization.dynamics
# #     Δt = 0.1 # seconds
# #     ϕ  = 2π*s.B0inho*Δt
    
# #     # Target Magnetization
# #     Mx_tar = M_tar[1,1]
# #     My_tar = M_tar[2,1]
# #     Mz_tar = M_tar[3,1]

# #     # Magnetization
# #     Mx = (m[2,end] - Mx_tar)^2/2
# #     My = (m[3,end] - My_tar)^2/2
# #     Mz = (m[4,end] - Mz_tar)^2/2

# #     c = (1 - (sin(ϕ)*Mx + cos(ϕ)*My + Mz))/s.Nspins

# #     return c
# # end
