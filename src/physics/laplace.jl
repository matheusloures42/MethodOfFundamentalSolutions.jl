struct DirichletType <: FieldType end

struct NeumannType <: FieldType end

struct LaplaceMedium{Dim,T} <: PhysicalMedium{Dim, 1} end

function greens(
    field::DirichletType, 
    medium::LaplaceMedium{2, T1},
    x::SVector{2, T2}, 
    outward_normal::SVector
) where {T1, T2}
    
    # Calculate the purely real scalar value
    val = -1 / (2π) * log(norm(x))
    
    # Return it wrapped in a 1x1 Static Matrix of type T
    return SMatrix{1, 1, T2}(val)
end

function greens_gradient(
    field::DirichletType, 
    medium::LaplaceMedium{2, T1}, 
    x::SVector{2, T2}, 
    outward_normal::SVector
) where {T1, T2}
    
    r2 = x[1]^2 + x[2]^2 + 1e-12 
    
    grad_x = -T2(1) / T2(2π) * (x[1] / r2)
    grad_y = -T2(1) / T2(2π) * (x[2] / r2)
    
   
    return SVector{2, T2}(grad_x, grad_y)
end


function laplace_M(medium::P, x_flat::AbstractVector{T1}, chi::AbstractVector{T2}) where {P<:LaplaceMedium, T1<:Real, T2<:Real}
    # CRUCIAL FIX: Automatically resolve whether to use Float64 or ForwardDiff.Dual
    NumType = promote_type(eltype(x_flat), eltype(chi))
    
    n_sensors = div(length(x_flat), 2)
    n_sources = div(length(chi), 2)
    
    # Preallocate using the dynamically resolved type
    M = zeros(NumType, n_sensors, n_sources)
    
    for i in 1:n_sensors
        # Use NumType instead of hardcoding Float64
        x_sens = SVector{2, NumType}(x_flat[2i-1], x_flat[2i])
        
        for j in 1:n_sources
            x_sour = SVector{2, NumType}(chi[2j-1], chi[2j])
            
            r_vec = x_sens - x_sour
            
            # This calls your greens signature, which handles NumType natively
            G_complex_mat = greens(DirichletType(), medium, r_vec, zero(r_vec))
            
            M[i, j] = real(G_complex_mat[1, 1])
        end
    end
    return M
end

"""
Analytical gradient tensor builder for Laplace.
Computes ∂M/∂x_sensor, where source locations are read from hyperparameter chi.
Returns a 3D array of shape (n_sensors, n_sources, 2)
"""

function laplace_grad_M(medium::P, x_flat::AbstractVector{T1}, chi::AbstractVector{T2}) where {P<:LaplaceMedium, T1<:Real, T2<:Real}
    # Apply the same type promotion logic here
    NumType = promote_type(eltype(x_flat), eltype(chi))
    
    n_sensors = div(length(x_flat), 2)
    n_sources = div(length(chi), 2)
    
    # Preallocate using NumType
    grad_M = zeros(NumType, n_sensors, n_sources, 2)
    
    for i in 1:n_sensors
        x_sens = SVector{2, NumType}(x_flat[2i-1], x_flat[2i])
        
        for j in 1:n_sources
            x_sour = SVector{2, NumType}(chi[2j-1], chi[2j])
            
            dx = x_sens[1] - x_sour[1]
            dy = x_sens[2] - x_sour[2]
            r2 = dx^2 + dy^2 + 1e-12 
            
            grad_M[i, j, 1] = -(1.0 / (2 * π)) * (dx / r2) 
            grad_M[i, j, 2] = -(1.0 / (2 * π)) * (dy / r2) 
        end
    end
    return grad_M
end