using ReinforcementLearning
using BlochSim

Base.@kwdef mutable struct BlochSimEnv <: AbstractEnv

    # The BlochSimEnv reward is the difference between samples A & B at time t
    reward::Float64 = 0.0
    M₀::Float64 = 0.0
    T₁::Float64 = 1000
    T₂::Float64 = 100
    Δf::Float64 = 3.75
    spin::Spin = Spin(M₀, T₁, T₂, Δf)
    elapsed_time::Float64 = 0.0  # ms
    max_t::Float64 = 10000  # ms
end

struct Pulse
    ω_x::Float64
    ω_y::Float64
end
RLBase.action_space(env::BlochSimEnv) = (Pulse, nothing)
RLBase.reward(env::BlochSimEnv) = env.reward
RLBase.state(env::BlochSimEnv) = (env.spin, env.reward)
RLBase.state_space(env::BlochSimEnv) = [false, true]
RLBase.is_terminated(env::BlochSimEnv) = env.elapsed_time <= env.max_t
RLBase.reset!(env::BlochSimEnv) = 
begin
    env.spin = Spin(env.M₀, env.T₁, env.T₂, env.Δf);
    env.reward = 0.0;
end

function (x::BlochSimEnv)(action)
    if isnothing(action)
        freeprecess!(env.spin, 0.1)
        # TODO: What is our reward?
        x.reward = 0.0
    elseif Pulse.typeof(action)
        # TODO: How are omegas related to the excitation?
        excite!(spin, InstantaneousRF(π/2))
        x.reward = 1.0
    else
        @error "unknown action of $action"
    end
end

env = BlochSimEnv()
RLBase.test_runnable!(env)
