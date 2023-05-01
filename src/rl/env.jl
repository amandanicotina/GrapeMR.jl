using ReinforcementLearning
using ReinforcementLearningBase
using BlochSim
using LinearAlgebra

Base.@kwdef mutable struct BlochSimEnv <: AbstractEnv

    # The BlochSimEnv reward is the difference between samples A & B at time t
    reward::Float64 = 0.0
    M₀::Float64 = 1.0  # Don't change this one
    T₁::Float64 = 1000
    T₂::Float64 = 100
    Δf::Float64 = 3.75
    spin::Spin = Spin(M₀, T₁, T₂, Δf)
    time_step::Float64 = 1  # ms per environment step
    elapsed_time::Float64 = 0.0  # ms
    max_t::Float64 = 10000  # ms
end

struct RfPulse
    ω_x::Float64
    ω_y::Float64
    Δt::Float64
end

RLBase.action_space(env::BlochSimEnv) = (RfPulse, nothing)
RLBase.reward(env::BlochSimEnv) = env.reward
RLBase.state(env::BlochSimEnv) = (env.spin.M)
RLBase.state_space(env::BlochSimEnv) = ArrayProductDomain([typemin(Float64) .. typemax(Float64), typemin(Float64) .. typemax(Float64), typemin(Float64) .. typemax(Float64)])
RLBase.is_terminated(env::BlochSimEnv) = env.elapsed_time <= env.max_t
RLBase.reset!(env::BlochSimEnv) = 
begin
    env.spin = Spin(env.M₀, env.T₁, env.T₂, env.Δf);
    env.reward = 0.0;
end

function (x::BlochSimEnv)(action)
    if isnothing(action)
        freeprecess!(env.spin, x.time_step)
        # TODO: What is our reward?
        x.reward = 0.0
    elseif isa(action, RfPulse)
        # TODO: How are omegas related to the excitation?
        excite!(spin, InstantaneousRF(π/2))
        freeprecess!(env.spin, x.time_step)
        x.reward = norm(normalize([0.5, 0.5, 0.5] - env.spin.M))
    else
        @error "unknown action of $action"
    end
    x.elapsed_time += x.time_step
end

env = BlochSimEnv()
RLBase.test_runnable!(env)
