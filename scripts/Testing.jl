using ChebyshevPicardIteration
using DifferentialEquations

# First Test
ϵ = 0.001
f(t,y) = cos(t .+ ϵ*y)

t0 = 0
tf = 256*π
y0 = 1

(ts,ys) = Integrate(f, y0, (t0, tf), 100, 200, 1e-14, 1000)

prob = ODEProblem((y,p,t) -> f(t,y), y0, (t0, tf))
sol = solve(prob, Tsit5(), reltol = 1e-8, abstol = 1e-8)

display(ys)
