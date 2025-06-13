################################################################################
#                            Abstract Pulse Shape Type                         #
################################################################################

"""
    PulseShape

Abstract type representing an RF pulse shape. Concrete subtypes implement 
different pulse generation strategies for control field design.
"""
abstract type PulseShape end

################################################################################
#                              Concrete Pulse Shapes                           #
################################################################################

"""
    Spline <: PulseShape

Cubic spline-based pulse shape using randomly sampled or user-defined control points.
Useful for generating smooth, interpolated RF envelopes.
"""
struct Spline <: PulseShape end

"""
    Hard <: PulseShape

Hard rectangular pulse shape with constant RF amplitude. Useful for quick excitation
or testing idealized spin behavior.
"""
struct Hard <: PulseShape end

"""
    Sinc <: PulseShape

Sinc-shaped pulse used for slice selection in MRI. Defined by flip angle and bandwidth.
B1x and B1y are phase-shifted by π/2.
"""
struct Sinc <: PulseShape end

"""
    Gaussian <: PulseShape

Pulse shape defined by a symmetric Gaussian profile. Often used for selective excitation 
or spectrally localized pulses.
"""
struct Gaussian <: PulseShape end

"""
    ZeroPulse <: PulseShape

A pulse shape with no RF excitation (B1x = B1y = 0). Useful for testing relaxation-only dynamics
or constructing periods of inactivity in a pulse sequence.
"""
struct ZeroPulse <: PulseShape end


"""
    BSSFP <: PulseShape

Balanced Steady-State Free Precession pulse type. Used in rapid imaging sequences with high SNR.
This shape assumes discrete pulses spaced over several TR periods.
"""
struct BSSFP <: PulseShape end


################################################################################
#                             Symbol-to-Trait Dispatch                         #
################################################################################

"""
    pulse_shape(::Val{Symbol}) -> PulseShape

Maps a `Symbol` to a concrete `PulseShape` type, enabling trait-based dispatching 
for RF pulse generation.

# Supported Symbols
- `:spline`   → [`Spline`](@ref)
- `:hard`     → [`Hard`](@ref)
- `:sinc`     → [`Sinc`](@ref)
- `:gaussian` → [`Gaussian`](@ref)
- `:bssfp`    → [`BSSFP`](@ref)
- `:zeropulse`    → [`ZeroPulse`](@ref)

# Example
```julia
julia> pulse_shape(Val(:spline))
Spline()

julia> pulse_shape(Val(:hard))
Hard()
```
"""
pulse_shape(::Val{:spline})    = Spline()
pulse_shape(::Val{:hard})      = Hard()
pulse_shape(::Val{:sinc})      = Sinc()
pulse_shape(::Val{:gaussian})  = Gaussian()
pulse_shape(::Val{:bssfp})     = BSSFP()
pulse_shape(::Val{:zeropulse}) = ZeroPulse()
