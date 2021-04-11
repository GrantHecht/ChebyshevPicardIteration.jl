
# Chebyshev Polynomial of the First Kind
function Tk(τk::AbstractFloat, k::Integer)
    return cos(k*acos(τk))
end

function T(τ::AbstractVector, N::Integer)
    # M + 1 (Number of nodes)
    Mp1 = length(τ)

    # Initialize Chebyshev Interpolation Matrix
    Tτ = SizedMatrix{N, Mp1}(ones(N, Mp1))

    # Fill Matrix
    Tτ[2, :] .= τ
    @inbounds Threads.@threads for col in 1:Mp1
        for row in 3:N
            Tτ[row,col] = Tk(τ[col], row - 1)
        end
    end
    return Tτ
end

function T!(Tτ::AbstractMatrix, τ::AbstractVector)
    dims = size(Tτ)

    # Fill Matrix
    Tτ[1, :] .= 1
    Tτ[2, :] .= τ
    @inbounds Threads.@threads for col in 1:dims[2]
        for row in 3:dims[1]
            Tτ[row,col] = Tk(τ[col], row - 1)
        end
    end
end

function ChebyshevInterpMat!(Cx::AbstractMatrix, τ::AbstractVector)
    dims = size(Cx)

    # Fill Matrix
    Cx[:,1] .= 1
    Cx[:,2] .= τ
    @inbounds Threads.@threads for col in 3:dims[2]
        for row in 1:dims[1]
            Cx[row,col] = Tk(τ[row], col - 1)
        end
    end
end
