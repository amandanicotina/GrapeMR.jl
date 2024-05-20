


function steady_state(s::GrapeMR.SteadyState)
    (TR, TE, T1, T2, M0) = (s.TR*1e3, s.TE*1e3, s.T1*1e3, s.T2*1e3, s.M_init[3])
    
    # Spin object
    spin = BlochSim.Spin(M0, T1, T2, s.B0inho)

    # Create RF pulse
    rf = InstantaneousRF(s.α)
    
    # Compute SS signal
    nTR = ceil(3*T1/TR)
    for i = 1:nTR
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
    #sig = BlochSim.signal(spin)

    return spin.M
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
    #sig = BlochSim.signal(spin)

    return spin.M
end


function steady_state_geometric(s::GrapeMR.SteadyState)
    (T1, T2, Δf, M0) = (s.T1, s.T2, s.B0inho, 1.0)
    (α, TR) = (s.α, s.TR)
    ϕ = 2π*Δf*TR # Offset angle
    β = 2*atan(tan(α/2) / cos(ϕ/2)) # β in rads

    Mxy = M0 / (cot(β/2) + (T1/T2)*tan(β/2))
    return Mxy
end

function steady_state_geometric_Mz(s::GrapeMR.SteadyState)
    (T1, T2, Δf, M0) = (s.T1, s.T2, s.B0inho, 1.0)
    (α, TR) = (s.α, s.TR)
    ϕ = 2π*Δf*TR # Offset angle
    β = 2*atan(tan(α/2) / cos(ϕ/2)) # β in rads

    Mz  = M0 * cos(ϕ/2) / (cot(β/2) + (T1/T2)*tan(β/2))
    return Mz
end

# Document functions like this:
"""
    save_grape_data(gp::GrapeMR.GrapeOutput; folder_path = pwd())

Save data related to Grape optimization into files organized in folders.

# Arguments
- `gp::GrapeMR.GrapeOutput`: Grape optimization output.
- `folder_path::String = pwd()`: Folder path where data will be saved.

# Example
```julia
save_grape_data(gp_output, folder_path="/path/to/folder")

If no path is provided, it saves the files inside the folder where the package was installed
folder name format : yyyy-mm-dd
"""