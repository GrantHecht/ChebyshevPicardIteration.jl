
function Integrate(f, x0, tspan, N, M, tol, imax)

    # Get length of state vector
    L = length(x0)

    # Transformation of independant variable
    ω₁ = 0.5*(tspan[2] + tspan[1])
    ω₂ = 0.5*(tspan[2] - tspan[1])

    g(τ,x) =   ω₂*f(ω₁ + ω₂*τ,x)

    # Compute Chebyshev-Gauss-Lobatto Nodes
    τs = SizedVector{M + 1}(zeros(M + 1))
    @inbounds for j in 0:M
        τs[j + 1] = -cos(j*π/M)
    end

    # Compute constant matricies
    Cx  = SizedMatrix{M + 1, N + 1}(zeros(M + 1, N + 1))
    ChebyshevInterpMat!(Cx, τs)
    CI1 = SizedMatrix{N + 1, N}(zeros(N + 1, N))
    initCI1!(CI1)
    Cf  = SizedMatrix{N, M + 1}(zeros(N, M + 1))
    initCf!(Cf, τs)
    C   = SizedMatrix{N + 1, M + 1}(zeros(N + 1, M + 1))
    mul!(C, CI1, Cf)

    # Initialize required matricies
    Xold  = SizedMatrix{M + 1, L}(zeros(M + 1, L))
    XoldT = Xold'
    Xnew  = SizedMatrix{M + 1, L}(zeros(M + 1, L))
    X0    = SizedMatrix{N + 1, L}(zeros(N + 1, L))
    X0[1,:] .= x0
    G     = SizedMatrix{M + 1, L}(zeros(M + 1, L))
    β     = SizedMatrix{N + 1, L}(zeros(N + 1, L))

    # Set Initial Guess for Xold (Could provide as input in future)
    Xold[1,:] .= x0

    # Begin integration loop
    stop = false
    eold = Inf
    enew = Inf
    iter = 0
    while !stop
        iter += 1
        println(iter)

        # Compute forcing function
        @inbounds for i in 1:M+1
            L == 1 ? G[i,1]  = g(τs[i], Xold[i]) :
                     G[i,:] .= g(τs[i], @view(XoldT[:, i]))
        end

        # Update Coefficients
        #β = X0 + C*G
        mul!(β,C,G)
        β .= β + X0

        # Update State
        #Xnew = Cx*β
        mul!(Xnew, Cx, β)

        # Compute error
        enew = 0
        for col in 1:L
            for row in 1:M+1
                temp = Xnew[row,col] - Xold[row,col]
                if temp > enew
                    enew = temp
                end
            end
        end

        if (enew <= tol && eold <= tol) || iter >= imax
            stop = true
        else
            Xold .= Xnew
            eold  = enew
        end
    end

    return (ω₁ .+ ω₂*τs, Xnew)
end

function initCI1!(CI1)
    # Get N
    dims = size(CI1)
    N = dims[2]

    # Initialize R and S
    R = Diagonal(SizedVector{N + 1}(zeros(N + 1)))
    S = SizedMatrix{N + 1, N}(zeros(N + 1, N))

    # Fill Matricies
    R[1,1] = 1
    R[2,2] = 0.5
    S[1,1] = 1
    S[1,2] = -0.25
    S[2,1] = 2
    S[2,3] = -1
    @inbounds for i in 3:N + 1
        R[i,i] = 1/(2*(i - 1))
        S[i,i - 1] = 1
        if i < N
            S[i,i + 1] = -1
        end
        if i < N + 1
            S[1,i] = (1/(2*(i - 2)) - 1/(2*i))*(-1)^i
        end
    end

    # Fill CI1
    #CI1 .= R*S
    lmul!(R,S)
    CI1 .= S
end

function initCf!(Cf, τs)
    # Get N and M
    dims = size(Cf)
    N = dims[1]
    M = dims[2] - 1

    # Set p variable depending on number of nodes and check Cf size for validity
    p = 1
    if M > N
        p = 2
    elseif N > M
        throw(ArgumentError("Improper size of Cf matrix!"))
    end

    # Initialize temporary matricies
    V  = Diagonal(SizedVector{N}(zeros(N)))
    W  = Diagonal(SizedVector{M + 1}(zeros(M + 1)))
    Tm = SizedMatrix{N, M + 1}(zeros(N, M + 1))

    # Fill T
    T!(Tm, τs)

    # Fill V
    @inbounds for i in 1:N
        # Get numerator
        i == 1 ? num = 1 : i == N ? num = p : num = 2

        # Set value
        V[i, i] = num/M
    end

    # Fill W
    @inbounds for i in 1:M + 1
        # Get val
        i == 1 ? val = 0.5 : i == M + 1 ? val = 0.5 : val = 1

        # Set value
        W[i, i] = val
    end

    # Compute Cf
    #Cf .= V*Tm*W
    rmul!(Tm,W)
    lmul!(V,Tm)
    Cf .= Tm
end
