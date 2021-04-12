using ChebyshevPicardIteration
using DifferentialEquations
using Plots

# First Test
ϵ = 0.001
f(t,y) = cos(t .+ ϵ*y)

t0 = 0
tf = 256*π
y0 = 1

# Integrate with Chebyshev-Picard Iteration
(ts,ys) = Integrate(f, y0, (t0, tf), 200, 200, 1e-20, 1000)

# Integrate with DifferentialEquations.jl
prob = ODEProblem((y,p,t) -> f(t,y), y0, (t0, tf))
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

# Analytic Solution
γ = ϵ/(1 + sqrt(1 - ϵ^2))
α = 2*ϵ/(1 - ϵ + sqrt(1 - ϵ^2))
β = tan(ϵ*y0/2)/(1 + α)
ϕ(t) = 0.5*(1 - γ*ϵ)*t
σ(t) = α*(sin(ϕ(t)) + β*cos(ϕ(t)))
y(t) = -γ*t + (2/ϵ)*atan((β + σ(t)*cos(ϕ(t)))/(1 + σ(t)*sin(ϕ(t))))

plot(sol.t,sol.u - map(y,sol.t))
plot!(ts,ys .- map(y,ts))

# Second Test
#g = 9.81
#k = 0.1
#m = 0.2
#f2(t,y) = [y[2], -k*y[1]/m + g]

#y0 = [5, 0]

# Integrate with Chebyshev-Picard Iteration
#(ts,ys) = Integrate(f2, y0, (t0, tf), 80, 80, 1e-14, 1000)

# Integrate with DifferentialEquations.jl
#prob = ODEProblem((y,p,t) -> f2(t,y), y0, (t0, tf))
#sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

#plot(sol)
#plot!(ts,ys)
