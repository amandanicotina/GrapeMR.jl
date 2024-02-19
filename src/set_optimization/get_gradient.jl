"""
    get_gradient(cf::ControlFields, iso::Magnetization, cost_function::String)

get_gradient
    # Input  
    - cf: (::InitialControlFields) - Control fields struct
    - iso: (::Magnetization) -
    - cost_function: (::String) - Key to the cost function dictionary

    # Output
    - ΔJx - 1xN matrix
"""
function get_gradient(cf::ControlFields, iso::Magnetization, H::Matrix, cost_function::String)
    χ     = backward_propagation(cf, iso, cost_function)
    t_arr = range(cf.t_control, 0.0, length = length(cf.B1x)+1)
    Δt    = diff(t_arr)[1]
    s     = iso.spin[1]
    M     = forward_propagation(cf, s)
    # Gradient 
    ΔJ = zeros(Float64, 1, length(cf.B1x))
    for i ∈ 1:length(cf.B1x)
        ΔJ[1,i] = transpose(χ[:,i+1])*H*Δt*M[:,i]
    end
    return ΔJ
end



"""
    update_control_field(cf::InitialControlFields, iso::Magnetization, cost_function::String, ϵ::Float64)

update_control_field
    # Input  
    - cf: (::InitialControlFields) - Control fields struct
    - iso: (::Magnetization) -
    - cost_function: (::String) - Key to the cost function dictionary
    - ϵ: (::Float64) - 

    # Output
    - Control Field - 1xN matrix
"""
function update_control_field(cf::ControlFields, iso::Magnetization, H::Matrix, cost_function::String, ϵ::Float64)
    ΔJ = get_gradient(cf, iso, H, cost_function)
    if H == Ix
        B = cf.B1x .+ ΔJ./ϵ
    else
        B = cf.B1y .+ ΔJ./ϵ
    end
    return B
end
