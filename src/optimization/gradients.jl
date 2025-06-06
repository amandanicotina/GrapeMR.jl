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
                H::AbstractMatrix{<:Real})
    @inbounds for i in 1:(size(M, 2) - 1)
        grad[1, i] = dot(
            transpose(view(χ, :, i + 1)),
            H,
            view(M, :, i + 1)
        )
    end
    return grad
end


"""
    update!(cf::ControlField, ∇xy::Tuple, ϵ::Float64)

Updates the control fields based on the calculated gradient and a learning rate.

# Arguments
- `cf::ControlField`: Control field struct to be updated.
- `∇xy::Tuple{Matrix{Float64}, Matrix{Float64}}`: Gradients for the x and y components of the field.
- `ϵ::Float64`: Learning rate for gradient descent.

# Returns
- `(u1x, u1y)`: Updated x and y control fields.
"""
function update!(cf::ControlField, ∇xy::Tuple{Matrix{Float64},Matrix{Float64}}, ϵ::Float64)
    u1x = cf.B1x .- ϵ .* ∇xy[1]
    u1y = cf.B1y .- ϵ .* ∇xy[2]
    return u1x, u1y
end