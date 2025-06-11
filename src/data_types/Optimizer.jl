"""
    AbstractOptimizer

Abstract supertype for all GRAPE optimization strategies.

Concrete subtypes include [`GradientDescent`](@ref) and [`BFGS`](@ref).
Used to dispatch to the appropriate optimization routine.
"""
abstract type AbstractOptimizer end

"""
    GradientDescent <: AbstractOptimizer

A GRAPE optimizer using first-order gradients and a polynomial decay
schedule for the step size `ϵ`.

Typically used in conjunction with [`GradientDescentConfig`](@ref), which
controls the decay parameters and number of iterations.

# Notes
- Simple and fast, but may converge slowly or get stuck in poor minima.
- Works well for most GRAPE applications with good initialization.
"""
struct GradientDescent <: AbstractOptimizer end

"""
    BFGS <: AbstractOptimizer

A GRAPE optimizer using the BFGS quasi-Newton method.

Stores an approximate inverse Hessian matrix to accelerate convergence
without computing second derivatives. Typically used with
[`BFGSConfig`](@ref).

# Notes
- Requires fewer iterations than gradient descent in many problems.
- More memory-intensive due to matrix storage.
- Not suitable if gradient noise is large or model is not smooth.
"""
struct BFGS <: AbstractOptimizer end

struct ManualGradientDescent <: AbstractOptimizer end

# -----------------------------------------------------------------------------

"""
    AbstractOptimizerConfig

Abstract supertype for all optimizer-specific configuration types used in the GRAPE algorithm.

Each optimizer (e.g., [`GradientDescent`](@ref), [`BFGS`](@ref)) requires its own configuration
with the appropriate hyperparameters. Subtypes of this abstract type define the structure
of those configurations.
"""
abstract type AbstractOptimizerConfig end

"""
    GradientDescentConfig <: AbstractOptimizerConfig

Configuration for the [`GradientDescent`](@ref) optimizer in GRAPE.

This configuration specifies a polynomial decay schedule for the step size `ϵ`, where
the learning rate decreases over iterations according to:
ϵ(k) = poly_start / (1 + k)^poly_degree

# Fields
- `poly_start::Float64`: Initial step size `ϵ₀`, controlling the starting update magnitude.
- `poly_degree::Int`: Degree of polynomial decay (how quickly `ϵ` decreases with iterations).
- `max_iter::Int`: Maximum number of GRAPE iterations.

# Example
```julia
config = GradientDescentConfig(0.75, 1, 300)
"""
struct GradientDescentConfig <: AbstractOptimizerConfig
    poly_start::Float64
    poly_degree::Int
    max_iter::Int
end

struct ManualGradientDescentConfig <: AbstractOptimizerConfig
    poly_start::Float64
    poly_degree::Int
    max_iter::Int
end

"""
    BFGSConfig <: AbstractOptimizerConfig

Configuration for the [`BFGS`](@ref) optimizer in GRAPE.

This configuration is minimal and includes only the number of iterations,
since BFGS uses an internal update rule to adapt step directions and magnitudes
based on an approximation of the inverse Hessian.

# Fields
- `max_iter::Int`: Maximum number of GRAPE iterations.

# Example
```julia
config = BFGSConfig(100)
"""
struct BFGSConfig <: AbstractOptimizerConfig
    max_iter::Int
end

Base.show(io::IO, grad_desc_congif::GradientDescentConfig) = 
    print(io, "GradientDescentConfig(poly_start = $(grad_desc_congif.poly_start), poly_degree = $(grad_desc_congif.poly_degree), max_iter = $(grad_desc_congif.max_iter))")
    
Base.show(io::IO, bfgs_congif::BFGSConfig) = 
    print(io, "BFGSConfig(max_iter = $(bfgs_congif.max_iter))")

