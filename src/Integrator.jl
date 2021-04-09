
function Integrate(f, x0, tspan, N, tol)

    # Get length of x0
    M = length(x0)

    # Transform independant variable
    ω1 = (tspan[2] + tspan[1])/2
    ω2 = (tspan[2] - tspan[1])/2

    # Define new function of tau
    g(τ,x) = ω2*f(ω1 + ω2*τ, x)'

    # Compute CGL nodes
    τs = SizedVector{N + 1}(zeros(N + 1))
    @inbounds for i in 0:N; τs[i] = -cos(i*π/N); end

    # Initialize Matricies and Tensors
    θx0  = SizedMatrix{N+1, M}(zeros(N+1, M))
    Cx   = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    Cα   = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    S    = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    R    = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    V    = SizedMatrix{N+1, N+1}(zeros(N+1, N+1))
    Xold = SizedMatrix{N+1, M}(zeros(N+1, M))
    Xnew = SizedMatrix{N+1, M}(zeros(N+1, M))
    G    = SizedMatrix{N+1, M}(zeros(N+1, M))
    β    = SizedMatrix{N+1, M}(zeros(N+1, M))

    # Fill Matricies
    R[1,1] = 1
    R[2,2] = 1/2
    S[1,1] = 1
    S[1,2] = -1/2
    S[2,1] = 1
    V[1,1] = 1/N
    V[N+1,N+1] = 1/N
    Cx[:,1] .= 1
    Cx[:,2] .= τs
    @inbounds for i in 2:N
        R[i+1, i+1] = 1/(2*i)
        V[i, i]     = 2/N
        S[i+1, i]   = 1
        @inbounds for j in 0:N
            Cx[j+1, i+1] = T(τs[j+1], i)
        end
        if i > 2
            S[1, i]   = (1/(i-2) - 1/i)*(-1)^i
            S[i-1, i] = -1
        end
    end
    mul!(Cα, Cx, V)
    mul!(V, S, Cα)
    mul!(Cα, R, V)

    # Begin Integration Loop
    e_old = 10*tol
    while e_old > tol

        # Compute Forcing Function
        @inbounds Threads.@threads for i in 0:N
            G[i,:] = g(τs[i+1], )
end
