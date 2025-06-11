"""
    gradient!(grad::AbstractMatrix, χ::Matrix{Float64}, M::Matrix{Float64}, H::AbstractMatrix)

In-place calculation of the gradient of the cost function with respect to the Hamiltonian.

# Arguments
- `grad::AbstractMatrix`: Preallocated 1×N matrix to store the gradient.
- `χ::Matrix{Float64}`: Adjoint state matrix.
- `M::Matrix{Float64}`: Forward magnetization matrix.
- `H::AbstractMatrix`: Hamiltonian matrix.

# Returns
- `grad`: Updated 1×N gradient matrix.
"""
function gradient!(grad::AbstractMatrix{<:Real},
                χ::AbstractMatrix{<:Real},
                M::AbstractMatrix{<:Real},
                H::AbstractMatrix{<:Real}
                )
    for i in 1:(size(M, 2) - 1)
        grad[1, i] = dot(
            transpose(view(χ, :, i + 1)),
            H,
            view(M, :, i + 1)
        )
    end
    return grad
end

"""
    gradient(χ::Matrix{Float64}, M::Matrix{Float64}, H::Matrix)

Calculates the gradient of the cost function with respect to the Hamiltonian for each time step.

# Arguments
- `χ::Matrix{Float64}`: Adjoint state matrix.
- `M::Matrix{Float64}`: Forward propagation matrix.
- `H::Matrix`: Hamiltonian matrix.

# Returns
- `grad::Matrix{Float64}`: Gradient of the cost function, as a 1xN matrix.
"""
function gradient(χ::Matrix{Float64},
            M::Matrix{Float64}, 
            H::AbstractMatrix{Int64}
            )
    grad = zeros(Float64, 1, size(M, 2) - 1)
    for i in 1:(size(M, 2) - 1)
        grad[1, i] = dot(
            transpose(view(χ, :, i + 1)),
            H,
            view(M, :, i + 1)
        )
    end
    return grad
end

"""
    update!(cf::AbstractControlField, ∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}, ϵ::Float64)

Computes the updated x and y components of the control field using gradient descent.

This function is compatible with both `ControlField` and `NormalizedControlField` types, allowing updates in either physical or normalized units.

# Arguments
- `cf::AbstractControlField`: The control field to be updated (can be normalized or physical).
- `∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}`: Tuple containing gradients of the cost function w.r.t. the B1x and B1y fields.
- `ϵ::Float64`: Learning rate used to scale the gradient step.

# Returns
- `(u1x, u1y)`: Tuple with updated B1x and B1y fields (not stored in `cf` in-place).
"""
function update!(cf::AbstractControlField, ∇xy::Tuple{Matrix{Float64},Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ .* ∇xy[1]
    u1y = cf.B1y .- ϵ .* ∇xy[2]
    return u1x, u1y
end
