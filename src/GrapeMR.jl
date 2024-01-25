module GrapeMR

using Plots
using BlochSim
using ExponentialUtilities

include("get_params.jl")
include("propagation/bloch_methods.jl")

export γ_¹H 
export InitialControlFields, Spins, OptimizationParams
export magnetization_ODE

end
