"""
    generate_control_field(::Spline; N, t_c, B1ref, rng)

Generates a cubic spline RF pulse.

# Keyword Arguments
- `N::Int`: Number of time points.
- `t_c::Float64`: Total control duration in seconds.
- `B1ref::Float64`: Reference RF amplitude.
- `rng`: Random number generator (optional, for reproducibility).

# Returns
- `ControlField`: A struct containing the spline-based control field.
"""
function generate_control_field(::Spline; N, t_c, B1ref, rng=Random.GLOBAL_RNG)
    len = 10
    spline_time = range(0.0, t_c; length=len)
    control_time = range(0.0, t_c; length=N)
    B1_vals = rand(rng, Float64, len)
    B1x = B1ref .* create_spline(spline_time, control_time, B1_vals)
    B1y = B1ref .* create_spline(spline_time, control_time, B1_vals)
    return build_control_field(B1x, B1y, B1ref, t_c)
end

"""
    generate_control_field(::Hard; N, t_c, B1ref)

Generates a rectangular (hard) RF pulse.
"""
function generate_control_field(::Hard; N, t_c, B1ref)
    B1x = B1ref .* ones(N)
    B1y = zeros(N)
    return build_control_field(B1x, B1y, B1ref, t_c)
end

"""
    generate_control_field(::Sinc; N, t_c, B1ref, α)

Generates a sinc-shaped RF pulse with phase offset.
"""
function generate_control_field(::Sinc; N, t_c, B1ref, α=π / 2)
    BW_Hz = 100.0
    t = range(0.0, t_c; length=N) .- t_c / 2
    flip = α / (2π * step(t))
    x = BW_Hz .* t
    B1x = (flip .* sinc.(x)) ./ 2π
    B1y = (flip .* sinc.(x .+ π / 2)) ./ 2π
    return build_control_field(B1x, B1y, B1ref, t_c)
end

"""
    generate_control_field(::Gaussian; N, t_c, B1ref)

Generates a Gaussian-shaped RF pulse.
"""
function generate_control_field(::Gaussian; N, t_c, B1ref)
    x = ((1:N) .- N / 2) ./ (N / 10)
    B1 = B1ref .* exp.(-0.5 .* x.^2)
    return build_control_field(B1, B1, B1ref, t_c)
end

"""
    generate_control_field(::BSSFP; N, t_c, B1ref, nTR, α, TR)

Generates a bSSFP-style pulse sequence. Assumes N is divisible by nTR.

# Errors
Throws if `N % nTR ≠ 0`.
"""
function generate_control_field(::BSSFP; N, t_c, B1ref, nTR, α, TR)
    if N % nTR != 0
        error("N must be divisible by nTR to evenly distribute time points across nTR.")
    end
    points_per_TR = N ÷ nTR
    Δt = TR / points_per_TR
    rf0 = (α / 2) / (2π * Δt)
    rf  = α / (2π * Δt)

    bSSFP_vec = zeros(N)
    for n in 1:nTR
        idx = (n - 1) * points_per_TR + 1
        bSSFP_vec[idx] = n == 1 ? rf0 : rf
    end

    return build_control_field(bSSFP_vec, zeros(N), B1ref, nTR * TR)
end

################################################################################
#                              Unified Interface                               #
################################################################################

"""
    generate_control_field(kind::Symbol = :spline; kwargs...)

Generates a control field given a pulse shape symbol. This is the main 
user-facing entry point for RF pulse generation.

# Supported `kind` values
- `:spline`, `:hard`, `:sinc`, `:gaussian`, `:bssfp`

# Keyword Arguments
Depends on the pulse shape. Common ones:
- `N`: Number of time steps (required)
- `t_c`: Control duration in seconds (required)
- `B1ref`: Reference amplitude of the RF field (required)

# Pulse-specific Keywords
- `:sinc`   → `α`: Flip angle in radians (default: π/2)
- `:bssfp`  → `nTR`: Number of TR repetitions (required)
              `α`: Flip angle in radians (required)
              `TR`: Repetition time in seconds (required)
- `:spline` → `rng`: Optional RNG seed for reproducibility

# Returns
- A `ControlField` struct

# Examples
```julia
# Spline-based RF control field:
generate_control_field(:spline; N=2000, t_c=0.5, B1ref=5.0)

# Hard (rectangular) pulse:
generate_control_field(:hard; N=2000, t_c=0.5, B1ref=5.0)

# Sinc-shaped pulse with a π/2 flip angle:
generate_control_field(:sinc; N=2000, t_c=0.5, B1ref=5.0, α=π/2)

# Gaussian-shaped pulse:
generate_control_field(:gaussian; N=2000, t_c=0.5, B1ref=5.0)

# Balanced SSFP pulse:
generate_control_field(:bssfp; N=2000, t_c=0.5, B1ref=5.0, nTR=10, α=π/2, TR=0.05)
```
"""
function generate_control_field(kind::Symbol = :spline; N::Int, t_c::Float64, B1ref::Float64, kwargs...)
    ps = pulse_shape(Val(kind))
    filtered = Dict(kwargs)

    if ps isa Spline
        rng = get(filtered, :rng, Random.GLOBAL_RNG)
        return generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref, rng=rng)

    elseif ps isa Hard
        return generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref)

    elseif ps isa Sinc
        α = get(filtered, :α, π / 2)
        return generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref, α=α)

    elseif ps isa Gaussian
        return generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref)

    elseif ps isa BSSFP
        nTR = get(filtered, :nTR, nothing)
        α   = get(filtered, :α, nothing)
        TR  = get(filtered, :TR, nothing)

        if any(isnothing, (nTR, α, TR))
            error("bSSFP requires keyword arguments: nTR, α, TR")
        end
        return generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref, nTR=nTR, α=α, TR=TR)

    else
        error("Unsupported pulse type: $kind")
    end
end
