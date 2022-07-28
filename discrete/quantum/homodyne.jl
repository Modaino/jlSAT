import QuantumOpticsBase: horner, hermite
import Polynomials: Polynomial

function compute_herms(fock::FockBasis)
    return Polynomial.(hermite.A(fock.N))
end

struct Homodyne
    herms::Vector{Polynomial{Float64}}
    coeffs::Vector{ComplexF64}
    function Homodyne(fock_a, fock_b)
        herms = compute_herms(fock_b)
        coeffs = zeros(ComplexF64, fock_a.N+1)
        return new(herms, coeffs)
    end
end

# Computes [f[1](x), ..., f[N](x)] where
# f[n](x) = C[n,k] HG[k](x) 
# HG[k](x) = <x|k>
# state = sum{n=1:N,k=1:K} C[n,k] |n>|k>
#       = sum{n=1:N} |n> int{dx} f[n](x) |x>
# Note: ptrace(state,1) = sum{n=1:N} int_{dx,dx'} f[n](x)f*[n](x') |x><x'|
function compute_coeffs!(hom::Homodyne, state::Ket, x::Number)
    hom.coeffs .= 0
    Gx = exp(-x^2/2) / π^(1/4)
    Hkx = n = k = 1
    for Cnk in state.data
        if n == 1
            Hkx = horner(hom.herms[k].coeffs, x)
        end
        hom.coeffs[n] += Cnk * Hkx * Gx
        n += 1
        if n > length(hom.coeffs)
            n = 1; k += 1
        end
    end
    return hom.coeffs
end

function homodyne_b(hom, state, rng, xmax=25)
    cdf = let ℋ=hom, ψ=state, w=rand(rng)
        x0 -> quadgk(x->norm(compute_coeffs!(ℋ,ψ,x))^2, -Inf, x0)[1] - w
    end
    return find_zero(cdf, (-xmax,xmax))
end

function homodyne_b!(hom, state, rng, xmax=25)
    x0 = homodyne_b(hom, state, rng, xmax)
    D = length(hom.coeffs)
    state.data[1 : D] .= compute_coeffs!(hom, state, x0)
    state.data[(D+1) : end] .= 0
    normalize!(state)
    return x0
end
