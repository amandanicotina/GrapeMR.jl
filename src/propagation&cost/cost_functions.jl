"""
Cost Functions 
introduce 1/Nspins
"""
function cost_function(iso::Isochromat, cf::Symbol)
    @match cf begin
        :euclidean_norm        => euclidean_norm(iso)
        :target_one_spin       => target_one_spin(iso)
        :target_steady_state   => target_steady_state(iso)
        :targetB0_steady_state => targetB0_steady_state(iso)
        :target_phase_encoding => target_phase_encoding(iso)
        :saturation_contrast              => saturation_contrast(iso)
        :saturation_contrast_Mx           => saturation_contrast(iso)
        :saturation_contrast_square       => saturation_contrast_square(iso)
        :saturation_contrast_steady_state => saturation_contrast_steady_state(iso)
        _                       => error("Cost function not defined")
    end
end

function calculate_cost_gradient(cost_func::Function, dict::Dict, vars::Vector{Num})
    dMx = Differential(vars[1])
    dMy = Differential(vars[2])
    dMz = Differential(vars[3])

    dCost_Mx = expand_derivatives(dMx(cost_func))
    dCost_My = expand_derivatives(dMy(cost_func))
    dCost_Mz = expand_derivatives(dMz(cost_func))
    Px = substitute.(dCost_Mx, (dict,))
    Py = substitute.(dCost_My, (dict,))
    Pz = substitute.(dCost_Mz, (dict,))

    cost_grad = [0.0, Px, Py, Pz]
    return cost_grad
 end
              
function euclidean_norm(iso::Isochromat)
    # Cost Function
    cost(Mx::Float64, My::Float64, Mz::Float64) = sqrt(Mx^2 + My^2 + Mz^2)

    # Numerical cost_values
    m = iso.magnetization.dynamics
    mx, my, mz = m[2,end], m[3,end], m[4,end]
    c = cost(mx, my, mz)

    # Calculating gradient
    vars = @variables Mx, My, Mz
    dict = Dict(vars .=> [mx, my, mz])
    c_grad = calculate_cost_gradient(cost, dict, vars)

    return c, c_grad
end

function target_one_spin(iso::Isochromat; M_tar = [0.0, 1.0, 0.0])

    # Target Magnetization
    Mx_tar = M_tar[1,1]
    My_tar = M_tar[2,1]
    Mz_tar = M_tar[3,1]

    # Magnetization
    mag = iso.magnetization.dynamics
    Mx = mag[2,end]
    My = mag[3,end]
    Mz = mag[4,end]

    c = ((Mx - Mx_tar)^2 + (My - My_tar)^2 + (Mz - Mz_tar)^2)/2

    return c
end

function target_phase_encoding(iso::Isochromat; M_tar = [0.0, 1.0, 0.0])
    s  = iso.spin
    m  = iso.magnetization.dynamics
    Δt = 0.1 # seconds
    ϕ  = 2π*s.B0inho*Δt
    
    # Target Magnetization
    Mx_tar = M_tar[1,1]
    My_tar = M_tar[2,1]
    Mz_tar = M_tar[3,1]

    # Magnetization
    Mx = (m[2,end] - Mx_tar)^2/2
    My = (m[3,end] - My_tar)^2/2
    Mz = (m[4,end] - Mz_tar)^2/2

    c = (1 - (sin(ϕ)*Mx + cos(ϕ)*My + Mz))/s.Nspins

    return c
end


function saturation_contrast(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        c = (1 - sqrt(m[4,end]*m[4,end] + 1e-15))/s.Nspins
    elseif s.target == "min"
        c = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)/s.Nspins
    end
    return c
end

function saturation_contrast_Mx(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        c = (1 - sqrt(m[2,end]*m[2,end] + 1e-15))/s.Nspins
    elseif s.target == "min"
        c = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)/s.Nspins
    end
    return c
end

function saturation_contrast_square(iso::Isochromat)
    c = 0.0;
    m = iso.magnetization.dynamics
    s = iso.spin

    if s.target == "max"
        c = (1 - sum(m[4,end]*m[4,end]))/s.Nspins
    elseif s.target == "min"
        c = sum(m[2:end,end].*m[2:end,end])/s.Nspins
    end
    return c
end

function target_steady_state(iso::Isochromat) 
    c = 0
    s = iso.spin
    if s.target == "min"
        # Steady State
        ss = steady_state_matrix(s)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
        # Magnetization
        mag = iso.magnetization.dynamics
        Mx  = mag[2,end]
        My  = mag[3,end]
        Mz  = mag[4,end]
        c = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)

    elseif s.target == "max"
        # Steady State
        ss = steady_state_matrix(s)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
        # Magnetization
        mag = iso.magnetization.dynamics
        Mx  = mag[2,end]
        My  = mag[3,end]
        Mz  = mag[4,end]
        c = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)

    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end

    return c/s.Nspins
end


function targetB0_steady_state(iso::Isochromat) 
    c = 0
    s = iso.spin

    # Steady State
    ss = steady_state_matrix(s)
    Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)

    # Magnetization
    mag = iso.magnetization.dynamics
    Mx  = mag[2,end]
    My  = mag[3,end]
    Mz  = mag[4,end]
    c = sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15)

    return c/s.Nspins
end

function saturation_contrast_steady_state(iso::Isochromat)
    c = 0
    s = iso.spin
    if s.target == "min"
        # Magnetization
        m = iso.magnetization.dynamics
        c = sqrt(sum(m[2:end,end].*m[2:end,end]) + 1e-15)/s.Nspins
    elseif s.target == "max"
        # Steady State
        ss = steady_state_matrix(s)
        Mx_ss, My_ss, Mz_ss = getproperty(ss, :x), getproperty(ss, :y), getproperty(ss, :z)
        # Magnetization
        m  = iso.magnetization.dynamics
        Mx = m[2,end]
        My = m[3,end]
        Mz = m[4,end]
        # c = (1 - sqrt((Mx - Mx_ss)^2 + (My - My_ss)^2 + (Mz - Mz_ss)^2 + 1e-15))/s.Nspins
        c = (1 - sqrt((Mz - Mz_ss)^2 + 1e-15))/s.Nspins
    else
        error(" $(s.target) is not a matching target. Valid targets are max or min")
    end

    return c
end



