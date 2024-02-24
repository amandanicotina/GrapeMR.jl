

"""
    get_gradient(cf::ControlFields, iso::Magnetization, cost_function::String)

get_gradient
    # Input  
    - cf:  (::InitialControlFields) - Control fields struct
    - iso: (::Magnetization) -
    - H:   (::Matrix) - Hamiltonian
    - cost_function: (::String) - Key to the cost function dictionary

    # Output
    - ΔJx - 1xN matrix
"""
function get_gradient(cf::ControlFields, iso::Magnetization, H::Matrix, cost_function::String)
    χ  = backward_propagation(cf, iso, cost_function)
    s  = iso.spin[1]
    M  = forward_propagation(cf, s)
    t  = range(cf.t_control, 0.0, length = length(cf.B1x)+1)
    Δt = diff(t)[1]
    # Gradient 
    ∇J = zeros(Float64, 1, length(cf.B1x))
    for i ∈ 1:length(cf.B1x)
        ∇J[1,i] = transpose(χ[:,i+1])*H*M[:,i]./2π
    end
    return ∇J
end



"""
    update_control_field(cf::InitialControlFields, iso::Magnetization, cost_function::String, ϵ::Float64)

update_control_field
    # Input  
    - cf:  (::InitialControlFields) - Control fields struct
    - iso: (::Magnetization) -
    - cost_function: (::String) - Key to the cost function dictionary
    - ϵ:   (::Float64) - 

    # Output
    - Control Field - 1xN matrix
"""
function update_control_field(cf::ControlFields, iso::Magnetization, H::Matrix, cost_function::String, ϵ::Float64)
    ∇J = get_gradient(cf, iso, H, cost_function)
    if H == Ix
        B = cf.B1x .- ϵ*∇J
    else
        B = cf.B1y .- ϵ*∇J
    end
    return B
end
