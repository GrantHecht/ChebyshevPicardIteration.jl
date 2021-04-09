using ChebyshevPicardIteration

# First Test
ϵ = 0.001
f(t,y) = cos(t .+ ϵ*y)

t0 = 0
tf = 256*π
y0 = 1

Integrate(f, y0, (t0, tf), 50, 1e-10, 1000)
