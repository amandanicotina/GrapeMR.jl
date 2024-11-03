
"""
    bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, Γ1::Float64, Γ2::Float64)

# Arguments
- B1x: (::Float64) - B1x step
- B1y: (::Float64) - B1x step
- Bz:  (::Float64) - B1x step
- Γ1:  (::Float64) - B1x step
- Γ2:  (::Float64) - B1x step

# Outputs
- Calculated 4x4 Bloch matrix
"""
function bloch_matrix(B1x::Float64, B1y::Float64, Bz::Float64, T1::Float64, T2::Float64)

    bloch_matrix = 
        SA[0.0    0.0   0.0    0.0;
         0.0   -1/T2  Bz    -B1y;
         0.0   -Bz   -1/T2   B1x;
         1/T1   B1y  -B1x   -1/T1] 
    
    return bloch_matrix
end

"""
    metal_bloch_matrix(B1x::Float32, B1y::Float32, Bz::Float32, Γ1::Float32, Γ2::Float32)

# Arguments
- B1x: (::Float32) - B1x step
- B1y: (::Float32) - B1x step
- Bz:  (::Float32) - B1x step
- Γ1:  (::Float32) - B1x step
- Γ2:  (::Float32) - B1x step
- Magnetization: (::MtlMatrix{Float32, Metal.SharedStorage}) - Magnetization matrix

# Outputs
- nothing

Note: Metal.jl kernels are not allowed to return values, so we need to pass the magnetization matrix as an argument
"""
function metal_bloch_matrix(B1x::Float32, B1y::Float32, Bz::Float32, T1::Float32, T2::Float32, Magnetization::MtlMatrix{Float32, Metal.SharedStorage})
    # We could also initialize a regular matrix and then feed it into the MtlArray constructor
    # TODO: Test this indirect usage of SA with regular matrix
    bloch_matrix = 
        SA[0.0f0    0.0f0   0.0f0    0.0f0;
         0.0f0   -1/T2  Bz    -B1y;
         0.0f0   -Bz   -1/T2   B1x;
         1/T1   B1y  -B1x   -1/T1]
    Magnetization .= bloch_matrix
    return
end


"""
forward_propagation

# Arguments  
- cf: (::ControlField) - Control fields struct
- s:  (::Spins) - Spin struct

# Outputs
- Magnetization vector 4xN
"""
function forward_propagation(cf::ControlField, s::Spins)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    M       = zeros(Float64, 4, length(cf.B1x)+1)
    M[:, 1] = [1.0, s.M_init[1], s.M_init[2], s.M_init[3]]
    
    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = cf.Bz .+ B0
    Bx = 2π*B1*cf.B1x
    By = 2π*B1*cf.B1y

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], s.T1,  s.T2)
        mul!(
            view(M, :, i+1), 
            exp(Δt*b_m),
            view(M, :, i)
        )
    end

    return M    
end

"""
metal_forward_propagation

Note: As we're using MtlArrays, we need to convert the Float64 to Float32
Additionally, all arguments to Metal kernels need to satisfy isbits

# Arguments
- cf: (::ControlField) - Control fields struct
- s:  (::Spins) - Spin struct

# Outputs
- nothing

Note: Metal.jl kernels are not allowed to return values, so we need to pass the magnetization matrix as an argument
"""
function metal_forward_propagation(cf_t_control::Float32, cf_B1x::Matrix{Float32}, cf_B1y::Matrix{Float32}, cf_Bz::Matrix{Float32}, s_B0inho::Float32, s_B1inho::Float32, s_T1::Float32, s_T2::Float32, s_M_init::Vector{Float32}, M::MtlMatrix{Float32, Metal.SharedStorage})
    Δt_arr  = MtlArray{Float32, 1, Metal.SharedStorage}(range(0.0, cf_t_control, length(cf_B1x)+1))
    M[:, 1] = MtlArray{Float32, 1, Metal.SharedStorage}([1.0, s_M_init[1], s_M_init[2], s_M_init[3]])
    
    B0 = Float32.(2π)*s_B0inho
    B1 = s_B1inho
    Bz = cf_Bz .+ B0
    Bx = Float32.(2π)*B1*cf_B1x
    By = Float32.(2π)*B1*cf_B1y

    for (i, Δt) ∈ enumerate(diff(Δt_arr))
        b_m = Metal.zeros(Float32, 4, 4; storage = Metal.SharedStorage)
        metal_bloch_matrix(Bx[i], By[i], Bz[i], s_T1,  s_T2, b_m)
        mul!(
            view(M, :, i+1), 
            exp(Δt*b_m),
            view(M, :, i)
        )
    end

    return
end

function test_forward_propagation(cf::ControlField, s::Spins)
    Δt_arr  = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt_diff = diff(Δt_arr)
    M       = zeros(Float64, 4, length(cf.B1x)+1)
    M[:, 1] = [1.0, s.M_init[1], s.M_init[2], s.M_init[3]];
    
    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = cf.Bz .+ B0
    Bx = 2π*B1.*cf.B1x
    By = 2π*B1.*cf.B1y 

    bloch_mat     = bloch_matrix.(Bx, By, Bz, s.T1, s.T2)
    exp_bloch_mat = [Matrix{Float64}(undef, 4, 4) for _ in 1:length(Δt_diff)]
    exp_bloch_mat = [exp(Δt * bloch_mat[i]) for (i, Δt) in enumerate(Δt_diff)]
    
    # Perform the forward propagation in a vectorized manner
    for i in 1:length(Δt_diff)
        M[:, i+1] = exp_bloch_mat[i] * M[:, i]
    end

    return M    
end


"""
backward_propagation

# Arguments  
- cf: (::ControlField) - Control fields struct
- iso: (::Isochromat) - Magnetization vector 4xN
- cost_function (::Function) - Cost Function gradient for adjoint state inital state

# Outputs
- Adjoint state 4xN
"""

function backward_propagation(cf::ControlField, iso::Isochromat, cost_grad::AbstractVector{Float64})
    t_arr      = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt         = diff(t_arr)
    back_steps = length(Δt)
    # TODO: refactor this as backward_propagation!(χ, cf, iso, cost_grad)
    χ          = zeros(Float64, 4, length(cf.B1x)+1)
    χ[:, end]  = cost_grad;
    s          = iso.spin

    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = 2π*cf.Bz .+ B0
    Bx = 2π*B1*cf.B1x
    By = 2π*B1*cf.B1y

    for i in back_steps:-1:1 
        b_m = bloch_matrix(Bx[i], By[i], Bz[i], s.T1, s.T2)
        bloch_matrix_adjoint = adjoint(b_m)
        mul!(
            view(χ, :, i),
            exp(Δt[i]*bloch_matrix_adjoint),
            view(χ, :, i+1)
        )
    end

    return round.(χ, digits = 5)
end

"""
backward_propagation

# Arguments  
- cf: (::ControlField) - Control fields struct
- iso: (::Isochromat) - Magnetization vector 4xN
- cost_function (::Function) - Cost Function gradient for adjoint state inital state

# Outputs
- nothing

Note: Metal.jl kernels are not allowed to return values, so we need to pass the adjoint magnetization matrix as an argument
"""

function metal_backward_propagation(cf_t_control::Float32, cf_B1x::Matrix{Float32}, cf_B1y::Matrix{Float32}, cf_Bz::Matrix{Float32}, s_B0inho::Float32, s_B1inho::Float32, s_T1::Float32, s_T2::Float32, χ::MtlMatrix{Float32, Metal.SharedStorage})
    t_arr      = range(0.0, cf_t_control, length(cf_B1x)+1)
    Δt         = Float32.(diff(t_arr))
    back_steps = length(Δt)
    # TODO: refactor this as backward_propagation!(χ, cf, iso, cost_grad)
    # χ          = Metal.zeros(Float32, 4, length(cf_B1x)+1; storage = Metal.SharedStorage)
    # χ[:, end]  = MtlArray{Float32, 1, Metal.SharedStorage}(cost_grad);

    B0 = Float32.(2π)*s_B0inho
    B1 = s_B1inho
    Bz = Float32.(2π)*cf_Bz .+ B0
    Bx = Float32.(2π)*B1*cf_B1x
    By = Float32.(2π)*B1*cf_B1y

    for i in back_steps:-1:1 
        b_m = Metal.zeros(Float32, 4, 4; storage = Metal.SharedStorage)
        metal_bloch_matrix(Bx[i], By[i], Bz[i], s_T1, s_T2, b_m)
        bloch_matrix_adjoint = adjoint(b_m)
        mul!(
            view(χ, :, i),
            exp(Δt[i]*bloch_matrix_adjoint),
            view(χ, :, i+1)
        )
    end
    # TODO: Skip for now, this fails compilation 
    # χ = round.(χ, digits = 5)
    return
end

function test_backward_propagation(cf::ControlField, iso::Isochromat, cost_grad::Vector{Float64})
    t_arr      = range(0.0, cf.t_control, length(cf.B1x)+1)
    Δt_diff    = diff(t_arr)
    back_steps = length(Δt_diff)
    χ          = zeros(Float64, 4, length(cf.B1x)+1)
    χ[:, end]  = cost_grad;
    s          = iso.spin

    B0 = 2π*s.B0inho
    B1 = s.B1inho
    Bz = 2π*cf.Bz .+ B0
    Bx = 2π*B1.*cf.B1x
    By = 2π*B1.*cf.B1y

    bloch_mat     = bloch_matrix.(Bx, By, Bz, s.T1, s.T2)
    bloch_mat_adj = adjoint.(bloch_mat)
    exp_bloch_mat = [Matrix{Float64}(undef, 4, 4) for _ in 1:back_steps]
    exp_bloch_mat = [exp(Δt_diff[i] * bloch_mat_adj[i]) for i in back_steps:-1:1]

    for i in back_steps:-1:1 
        χ[:, i] = exp_bloch_mat[i]*χ[:, i+1]
    end

    return round.(χ, digits = 5)
end
