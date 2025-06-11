"""
    normalize_rf!(B1x::Vector{Float64}, B1y::Vector{Float64}, B1ref::Float64) -> (B1x, B1y)

Normalizes the RF field vectors `B1x` and `B1y` such that the maximum magnitude 
of the complex-valued RF amplitude √(B1x² + B1y²) is equal to `B1ref`.

The function modifies the input vectors in-place.

# Arguments
- `B1x::Vector{Float64}`: x-component of the RF pulse in Hz.
- `B1y::Vector{Float64}`: y-component of the RF pulse in Hz.
- `B1ref::Float64`: Reference RF amplitude (in Hz) used to scale the pulse.

# Returns
- `(B1x, B1y)`: The normalized RF vectors (same objects, modified in-place).

# Notes
- If the maximum magnitude of the original RF is zero, the pulse is left unchanged.
- This ensures that the peak amplitude of the resulting RF pulse is `B1ref`.

"""
function normalize_rf!(B1x::Vector{Float64}, B1y::Vector{Float64}, B1ref::Float64)
    mag = sqrt.(B1x.^2 .+ B1y.^2)
    max_mag = maximum(mag)
    scale = max_mag > 0 ? max_mag : 1.0
    B1x .= (B1x ./ scale) .* B1ref
    B1y .= (B1y ./ scale) .* B1ref
    return B1x, B1y
end

"""
    create_spline(spline_time, control_time_vals, B1_vals; rng=Random.GLOBAL_RNG)

Creates a cubic spline interpolation from the specified control points and
evaluates it at given control time values.

# Arguments
- `spline_time::AbstractVector`: Time vector for knot positions.
- `control_time_vals::AbstractVector`: Time vector to evaluate the spline on.
- `B1_vals::AbstractVector`: Amplitude values for interpolation knots.

# Keywords
- `rng`: Random number generator (default: `Random.GLOBAL_RNG`).

# Returns
- Vector of interpolated values at `control_time_vals`.
"""
function create_spline(
    spline_time::AbstractVector,
    control_time_vals::AbstractVector,
    B1_vals::AbstractVector;
    rng = Random.GLOBAL_RNG,
)
    spline = CubicSpline(spline_time, B1_vals)
    return map(t -> spline(t), control_time_vals)
end


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

    B1_vals_x = rand(rng, len)
    B1_vals_y = rand(rng, len)

    B1x_raw = create_spline(spline_time, control_time, B1_vals_x)
    B1y_raw = create_spline(spline_time, control_time, B1_vals_y)

    normalize_rf!(B1x_raw, B1y_raw, B1ref)
    return build_control_field(B1x_raw, B1y_raw, B1ref, t_c)
end

"""
    generate_control_field(::Hard; N, t_c, B1ref)

Generates a rectangular (hard) RF pulse.
"""
function generate_control_field(::Hard; N, t_c, B1ref)
    B1x = fill(B1ref, N)
    B1y = zeros(N)
    # normalize_rf!(B1x, B1y, B1ref)
    return build_control_field(B1x, B1y, B1ref, t_c)
end

"""
    generate_control_field(::Sinc; N, t_c, B1ref, α)

Generates a sinc-shaped RF pulse with optional phase offset.
"""
function generate_control_field(::Sinc; N, t_c, B1ref, α=π/2)
    BW_Hz = 100.0
    t = range(0.0, t_c; length=N) .- t_c / 2
    flip = α / (2π * step(t))
    x = BW_Hz .* t
    B1x = (flip .* sinc.(x)) ./ 2π
    B1y = (flip .* sinc.(x .+ π / 2)) ./ 2π
    normalize_rf!(B1x, B1y, B1ref)
    return build_control_field(B1x, B1y, B1ref, t_c)
end

"""
    generate_control_field(::Gaussian; N, t_c, B1ref)

Generates a Gaussian-shaped RF pulse.
"""
function generate_control_field(::Gaussian; N, t_c, B1ref)
    x = ((1:N) .- N / 2) ./ (N / 10)
    B1 = exp.(-0.5 .* x.^2)
    B1x = copy(B1)
    B1y = copy(B1)
    normalize_rf!(B1x, B1y, B1ref)
    return build_control_field(B1x, B1y, B1ref, t_c)
end

"""
    ZeroPulse <: PulseShape

A pulse shape with zero RF amplitude (`B1x = B1y = 0`) across its duration.

Useful for simulating relaxation or free precession without RF excitation, or for inserting delay periods in pulse sequences.

# Example
```julia
cf = generate_control_field(ZeroPulse(); N, t_c, B1ref)
"""
function generate_control_field(::ZeroPulse;  N, t_c, B1ref)
    B1x = zeros(N)
    B1y = zeros(N)
    Bz = zeros(N)
    return build_control_field(B1x, B1y, B1ref, t_c)
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

    B1x = zeros(N)
    for n in 1:nTR
        idx = (n - 1) * points_per_TR + 1
        B1x[idx] = n == 1 ? rf0 : rf
    end
    B1y = zeros(N)
    normalize_rf!(B1x, B1y, B1ref)
    return build_control_field(B1x, B1y, B1ref, nTR * TR)
end

"""
    generate_control_field(kind::Symbol = :spline; kwargs...) -> AbstractControlField

Generates a control field (RF pulse) based on the selected pulse shape.

By default, the generated pulse is normalized to **unitless amplitude and time** using `B1ref`
as the sole normalization reference. This yields a `NormalizedControlField` suitable for 
optimization routines like GRAPE, where time is scaled as `t ⋅ B1ref`.

# Supported `kind` values
- `:spline`: Cubic spline pulse with randomized control points.
- `:hard`: Rectangular pulse (constant amplitude).
- `:sinc`: Sinc-shaped pulse with optional phase offset.
- `:gaussian`: Gaussian-shaped pulse.
- `:bssfp`: Balanced SSFP-style pulse train.

# Common Keyword Arguments
- `t_c::Float64`: Total control duration in seconds (required)
- `B1ref::Float64`: Reference RF amplitude in Tesla (required)
- `Δt_target::Float64`: Approximate target time resolution in seconds (default: `1e-4`)
- `normalize::Bool`: Whether to return a normalized control field (default: `true`)

# Pulse-Specific Keywords
- `:spline` → `rng`: Random number generator for reproducibility
- `:sinc` → `α`: Flip angle in radians (default: `π/2`)
- `:bssfp` → 
    - `nTR`: Number of TR repetitions (required)  
    - `α`: Flip angle in radians (required)  
    - `TR`: Repetition time in seconds (required)

# Returns
- `NormalizedControlField` (if `normalize=true`)  
  with time normalized as `t ⋅ B1ref` and amplitude as `B1 / B1ref`
- `ControlField` (if `normalize=false`)  
  in physical SI units (Tesla, seconds)

# Examples
```julia
cf = generate_control_field(:spline; t_c=0.5, B1ref=5.0)
cf_phys = generate_control_field(:hard; t_c=0.5, B1ref=5.0, normalize=false)
"""
function generate_control_field(kind::Symbol = :spline; 
    t_c::Float64,
    Δt_target::Float64=1e-4, 
    B1ref::Float64, 
    normalize::Bool=true,
    kwargs...
)
    N = max(10, ceil(Int, t_c / Δt_target))
    @debug "Auto-selected N = $N for t_c = $t_c and Δt_target = $Δt_target"

    ps = pulse_shape(Val(kind))
    filtered = Dict(kwargs)

    cf = if ps isa Spline
        rng = get(filtered, :rng, Random.GLOBAL_RNG)
        generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref, rng=rng)

    elseif ps isa Hard
        generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref)

    elseif ps isa Sinc
        α = get(filtered, :α, π / 2)
        generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref, α=α)

    elseif ps isa Gaussian
        generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref)

    elseif ps isa ZeroPulse
        generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref)

    elseif ps isa BSSFP
        nTR = get(filtered, :nTR, nothing)
        α   = get(filtered, :α, nothing)
        TR  = get(filtered, :TR, nothing)

        if any(isnothing, (nTR, α, TR))
            error("bSSFP requires keyword arguments: nTR, α, TR")
        end
        generate_control_field(ps; N=N, t_c=t_c, B1ref=B1ref, nTR=nTR, α=α, TR=TR)

    else
        error("Unsupported pulse type: $kind")
    end

    return normalize ? normalize_control_field(cf) : cf
end

