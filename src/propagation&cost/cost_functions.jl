

############################################################################################
#                                     Cost Functions                                       #
############################################################################################

norm_cost(x::AbstractVector{T}, Nspins) where {T} = sqrt(sum(abs2, x) + 1e-12) / Nspins
max_cost(x::AbstractVector{T}, Nspins) where {T} = (1 - sqrt(sum(abs2, x) + 1e-12)) / Nspins

function euclidean_norm(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    x = SVector(
        m[2, end], 
        m[3, end], 
        m[4, end]
    )
    val = norm_cost(x, s.Nspins)
    grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    return val, vcat(zero(eltype(grad)), grad)
end

function spin_target(iso::Isochromat; target::AbstractVector = [0.0, 1.0, 0.0])
    s = iso.spin
    m = iso.magnetization.dynamics
    x = SVector(
        (m[2, end] - target[1]),
        (m[3, end] - target[2]),
        (m[4, end] - target[3])
    )
    val = norm_cost(x, s.Nspins)
    grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    return val, vcat(zero(eltype(grad)), grad)
end

function saturation_contrast(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    if s.target == "max"
        x = SVector(
            zero(eltype(m[4, end])), 
            zero(eltype(m[4, end])), 
            m[4, end]
        )
        val = max_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(max_cost, s.Nspins), x)
        
    elseif s.target == "min"
        x = SVector(
            m[2, end], 
            m[3, end], 
            m[4, end]
        )
        val = norm_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return val, vcat(zero(eltype(grad)), grad)
end

function saturation_contrast_Mx(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    if s.target == "max"
        x = SVector(
            m[2, end], 
            zero(eltype(m[2, end])), 
            zero(eltype(m[2, end])), 
        )
        val = max_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(max_cost, s.Nspins), x)
    elseif s.target == "min"
        x = SVector(
            m[2, end], 
            m[3, end], 
            zero(eltype(m[2, end])) # It does not force the z-component to be zero
        )
        val = norm_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return val, vcat(zero(eltype(grad)), grad)
end

function saturation_contrast_Mtrans(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    if s.target == "max"
        x = SVector(
            m[2, end], 
            m[3, end], 
            zero(eltype(m[2, end]))
        )
        val = max_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(max_cost, s.Nspins), x)
    elseif s.target == "min"
        x = SVector(
            m[2, end], 
            m[3, end], 
            zero(eltype(m[2, end]))
        )
        val = norm_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return val, vcat(zero(eltype(grad)), grad)
end

function target_steady_state(iso::Isochromat) 
    s = iso.spin
    m = iso.magnetization.dynamics
    ss = steady_state_matrix(iso)
    Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

    x = SVector(
        (m[2, end] - Mx_ss),
        (m[3, end] - My_ss),
        (m[4, end] - Mz_ss)
    )

    val = norm_cost(x, s.Nspins)
    grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    return val, vcat(zero(eltype(grad)), grad)
end

function saturation_contrast_steady_state(iso::Isochromat)
    s = iso.spin
    m = iso.magnetization.dynamics
    if s.target == "max"
        ss = steady_state_matrix(iso)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
        x = SVector(
            (m[2, end] - Mx_ss),
            (m[3, end] - My_ss),
            (m[4, end] - Mz_ss)
        )
    
        val = norm_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)

    elseif s.target == "min"
        x = SVector(
            m[2, end], 
            m[3, end], 
            m[4, end]
        )
        val = norm_cost(x, s.Nspins)
        grad = ForwardDiff.gradient(Base.Fix2(norm_cost, s.Nspins), x)
    else
        error("Invalid target $(s.target). Valid targets are 'max' or 'min'.")
    end

    return val, vcat(zero(eltype(grad)), grad)
end

# function cos_func_do_dawid(iso::Isochromat)
#     vars = @variables Mx, My, Mz
#     cost = (sin(Mx) - cos(My))*atan(Mz)
#     return get_cost_and_gradients(iso, cost, vars)
# end


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


