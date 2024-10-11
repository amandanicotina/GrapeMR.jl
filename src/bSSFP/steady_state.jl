"""
    steady_state(s::GrapeMR.SteadyState)

Calculates the steady-state signal of a bSSFP sequence for a given set of parameters. The function calculates the magnetization evolution 
using repeated RF excitations and free precession steps until the steady state is reached.

    # Arguments
    - `s::GrapeMR.SteadyState`: A struct containing TR (Repetition Time), TE (Echo Time), T1 (Spin-lattice relaxation time), 
                                T2 (Spin-spin relaxation time), α (Flip angle in radians), Δϕ (Phase cycling increment), 
                                M_init (Initial magnetization vector), and B0inho (B0 inhomogeneity).
    # Returns
    - `signal::Vector{ComplexF64}`: The complex steady-state signal of the spin system.
"""
function steady_state(s::GrapeMR.SteadyState)
    (TR, TE, T1, T2, M0) = (s.TR*1e3, s.TE*1e3, s.T1*1e3, s.T2*1e3, s.M_init[3])
    
    # Spin object
    spin = BlochSim.Spin(M0, T1, T2, s.B0inho)

    # Create RF pulse
    rf = InstantaneousRF(s.α)
    
    # Compute SS signal
    nTR = ceil(Int, 3*T1/TR)
    for _ in 1:nTR
        excite!(spin, rf)
        freeprecess!(spin, TR)
        rf = InstantaneousRF(s.α, rf.θ + s.Δϕ)
    end
    excite!(spin, rf)
    freeprecess!(spin, TE)
    sig = BlochSim.signal(spin) * cis(-rf.θ)

    return sig
end


function steady_state_matrix(s::GrapeMR.SteadyState)
    # Create Spin object
    (TR, TE, T1, T2, M0) = (s.TR*1e3, s.TE*1e3, s.T1*1e3, s.T2*1e3, s.M_init[3])
    I  = BlochSim.I
    # Spin objects
    spin = BlochSim.Spin(M0, T1, T2, s.B0inho)
    spin_phase_cycle = BlochSim.Spin(M0, T1, T2, s.B0inho - (s.Δϕ / 2π / (TR/1000)))

    # Create RF pulse
    rf = InstantaneousRF(s.α)

    # Compute SS signal
    (R,)   = excite(spin_phase_cycle, rf)
    (A, B) = freeprecess(spin_phase_cycle, TR)
    M      = (I - R * A) \ (R * B)

    copyto!(spin.M, M)
    freeprecess!(spin, TE)
    sig = BlochSim.signal(spin)

    return sig
end


function steady_state_matrix(Tr, Te, T₁, T₂, M₀, B0inho, α, Δϕ)
    # Create Spin object
    (TR, TE, T1, T2, M0, B0) = (Tr*1e3, Te*1e3, T₁*1e3, T₂*1e3, M₀, B0inho[])
    I  = BlochSim.I
    # Spin objects
    spin = BlochSim.Spin(M0, T1, T2, B0)
    spin_phase_cycle = BlochSim.Spin(M0, T1, T2, B0 - (Δϕ / 2π / (TR/1000)))

    # Create RF pulse
    rf = InstantaneousRF(α)

    # Compute SS signal
    (R,)   = excite(spin_phase_cycle, rf)
    (A, B) = freeprecess(spin_phase_cycle, TR)
    M      = (I - R * A) \ (R * B)

    copyto!(spin.M, M)
    freeprecess!(spin, TE)

    return spin.M
end

function steady_state_matrix(iso::Isochromat)
    s = iso.spin
    (TR, TE, T1, T2, M0) = (s.TR*1e3, s.TE*1e3, s.T1*1e3, s.T2*1e3, s.M_init[3])
    I  = BlochSim.I
    # Spin objects
    spin = BlochSim.Spin(M0, T1, T2, s.B0inho)
    spin_phase_cycle = BlochSim.Spin(M0, T1, T2, s.B0inho - (s.Δϕ / 2π / (TR/1000)))

    # Create RF pulse
    rf = InstantaneousRF(s.α)

    # Compute SS signal
    (R,)   = excite(spin_phase_cycle, rf)
    (A, B) = freeprecess(spin_phase_cycle, TR)
    M      = (I - R * A) \ (R * B)

    copyto!(spin.M, M)
    freeprecess!(spin, TE)

    return spin.M
end

"""
    steady_state_geometric(s::GrapeMR.SteadyState)

Calculates the transverse steady-state magnetization using a geometric solution for a bSSFP sequence. Uses a geometric approach to 
derive the transverse component of the steady-state magnetization.

    # Arguments
    - `s::GrapeMR.SteadyState`: Struct containing the sequence parameters and spin system properties.

    # Returns
    - `Mxy::Vector{Float64}`: The steady-state transverse magnetization magnitude `Mxy`
"""
function steady_state_geometric(s::GrapeMR.SteadyState)
    (T1, T2, Δf, M0) = (s.T1, s.T2, s.B0inho, 1.0)
    (α, TR, Δϕ) = (s.α, s.TR, s.Δϕ)
    ϕ = 2π*Δf*TR # Offset angle
    β = 2*atan(tan(α/2) / cos(ϕ/2)) # β in rads

    Mxy = M0 / (cot(β/2) + (T1/T2)*tan(β/2))
    return Mxy
end

"""
    steady_state_geometric_Mz(s::GrapeMR.SteadyState)

Calculates the longitudinal steady-state magnetization using a geometric solution for a bSSFP sequence. Uses a geometric approach to 
derive the longitudinal component of the steady-state magnetization.

    # Arguments
    - `s::GrapeMR.SteadyState`: Struct containing the sequence parameters and spin system properties.

    # Returns
    - `Mz::Vector{Float64}`: The steady-state longitudinal magnetization magnitude `Mz`
"""
function steady_state_geometric_Mz(s::GrapeMR.SteadyState)
    (T1, T2, Δf, M0) = (s.T1, s.T2, s.B0inho, 1.0)
    (α, TR, Δϕ) = (s.α, s.TR, s.Δϕ) 
    ϕ = 2π*Δf*TR # Offset angle
    β = 2*atan(tan(α/2) / cos(ϕ/2)) # β in rads

    Mz  = (M0 * cot(β/2)) / (cot(β/2) + (T1/T2)*tan(β/2))
    return Mz
end



