
# Chebyshev Polynomial of the First Kind
function T(τ::AbstractFloat, k::Integer)
    return cos(k*acos(τ))
end
