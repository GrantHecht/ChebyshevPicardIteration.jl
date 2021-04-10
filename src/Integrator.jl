
function Integrate(f, x0, tspan, N, tol, imax)

    # Get length of x0
    M = length(x0)

    # Transform independant variable
    ω₁ = (tspan[2] + tspan[1])/2
    ω₂ = (tspan[2] - tspan[1])/2

    # Define new function of tau
    g(τ,x) = ω₂*f(ω₁ + ω₂*τ, x)

    # Compute CGL nodes
    τs = SizedVector{N + 1}(zeros(N + 1))
    @inbounds for i in 0:N; τs[i + 1] = -cos(i*π/N); end

    # Initialize Constant Matricies
    Cx = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    Cα = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    X0 = SizedMatrix{N+1, M}(zeros(N+1, M))
    initCx!(Cx, τs)
    initCα!(Cα, Cx)
    X0[1,:] .= x0

    # Initialize Requred Matrices
    Xold = SizedMatrix{N+1, M}(zeros(N+1, M))
    Xnew = SizedMatrix{N+1, M}(zeros(N+1, M))
    G    = SizedMatrix{N+1, M}(zeros(N+1, M))
    β    = SizedMatrix{N+1, M}(zeros(N+1, M))

    # Create Transposed Pointer to Xold
    XoldT = Xold'

    # Begin Integration Loop
    Xold[1,:] .= x0
    stop  = false
    e_old = Inf
    e_new = Inf
    i = 0
    while !stop

        # Increment Iteration Counter
        i += 1

        # Compute Forcing Function
        @inbounds for i in 0:N
            M == 1 ? G[i+1,1]  = g(τs[i+1], Xold[i+1]) :
                     G[i+1,:] .= g(τs[i+1], XoldT[:, i+1])
        end

        # Update Coefficients
        mul!(β, Cα, G)
        β .+= X0

        # Update State
        mul!(Xnew, Cx, β)

        # Compute Error
        e_new = 0
        for i in 0:N
            for j in 1:M
                temp =  abs(Xnew[i+1,j] - Xold[i+1,j])
                if temp > e_new; e_new = temp; end
            end
        end
        e_new = sqrt(e_new)/((N+1)*M)

        # Check for Convergence
        if (e_new <= tol && e_old <= tol) || i >= imax
            stop = true
        else
            Xold .= Xnew
            e_old = e_new
        end
    end

    return (ω₁ .+ ω₂.*τs, Xnew)
end

# Initialize Chebyshev Polynomial Matrix Cx
function initCx!(Cx::AbstractArray, τs::AbstractArray)
    dims = size(Cx)
    @inbounds for j in 1:dims[2]
        if j > 2
            for i in 1:dims[1]
                Cx[i,j] = T(τs[i], j - 1)
            end
        elseif j == 2
            for i in 1:dims[1]
                Cx[i,j] = τs[i]
            end
        else
            for i in 1:dims[1]
                Cx[i,j] = 1
            end
        end
    end
end

# Initialize Chebyshev Matrix Cα
function initCα!(Cα::AbstractArray, Cx::AbstractArray)
    dims = size(Cα)
    N = dims[2]

    # Initialize Temporary Matricies
    S = SizedMatrix{N, N}(zeros(N, N))
    R = SizedMatrix{N, N}(zeros(N, N))
    V = SizedMatrix{N, N}(zeros(N, N))

    # Fill Matricies
    N -= 1
    R[1,1]     = 1
    R[2,2]     = 1/2
    S[1,1]     = 1
    S[1,2]     = -1/2
    S[2,1]     = 1
    V[1,1]     = 1/N
    V[N+1,N+1] = 1/N
    @inbounds Threads.@threads for i in 2:N
        R[i+1, i+1] = 1/(2*i)
        V[i, i]     = 2/N
        S[i+1, i]   = 1
        if i > 2
            S[1, i]   = (1/(i-2) - 1/i)*(-1)^i
            S[i-1, i] = -1
        end
    end
    mul!(Cα, Cx, V)
    mul!(V, S, Cα)
    mul!(Cα, R, V)
end
