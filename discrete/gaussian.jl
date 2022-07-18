import Random: AbstractRNG, randn
import LinearAlgebra: I, mul!, BLAS.gemm!, BLAS.ger!

struct State
    μ::Vector{Float64}
    Σ::Matrix{Float64}
    cache::Vector{Float64}
    State(μ::Vector{<:Real}, Σ::Matrix{<:Real}) = new(μ, Σ, zero(μ))
end

Vacuum(N::Int) = State(zeros(2N), Matrix(0.5I, 2N, 2N))

function vacuum_project!(Γ::State, m::Int)
    Γ.μ[2m-1] = Γ.μ[2m] = 0
    for i in eachindex(Γ.μ)
        Γ.Σ[2m-1,i] = Γ.Σ[2m,i] = Γ.Σ[i,2m-1] = Γ.Σ[i,2m] = 0
    end
    Γ.Σ[2m-1,2m-1] = Γ.Σ[2m,2m] = 0.5
    return Γ
end

function squeeze_project!(Γ::State, m::Int, λ::Real)
    vacuum_project!(Γ, m)
    Γ.Σ[2m-1,2m-1] *= λ
    Γ.Σ[2m,2m] /= λ
end

struct Symplectic
    S::Matrix{Float64}
    cache::Matrix{Float64}
    Symplectic(S::Matrix{<:Real}) = new(S, zero(S))
end
Beamsplitter(R::Real) = let r=sqrt(R), t=sqrt(1-R)
    return Symplectic([
        t  0 -r  0 ;
        0  t  0 -r ;
        r  0  t  0 ;
        0  r  0  t ])
end

function update!(T::Symplectic, Γ::State)
    mul!(Γ.cache, T.S, Γ.μ)
    Γ.μ .= Γ.cache
    gemm!('N', 'T', 1.0, Γ.Σ, T.S, 0.0, T.cache)
    mul!(Γ.Σ, T.S, T.cache)
    return Γ
end

function displace_q!(Γ::State, m::Int, q::Real)
    Γ.μ[2m-1] += q
    return Γ
end

function homodyne_q!(Γ::State, m::Int, rng::AbstractRNG)
    ind = 2m-1
    q = Γ.μ[ind]
    σ² = Γ.Σ[ind,ind]
    z = randn(rng)
    for i in eachindex(Γ.μ)
        Γ.cache[i] = Γ.Σ[i,ind] / sqrt(σ²)
        Γ.μ[i] += z * Γ.cache[i]
    end
    ger!(-1.0, Γ.cache, Γ.cache, Γ.Σ)
    vacuum_project!(Γ, m)
    return sqrt(σ²) * z + q
end

function homodyne_q!(Γ::State, m::Int, rng::ThermalNoise)
    q = Γ.μ[2m-1]
    vacuum_project!(Γ, m)
    return q + randn(rng)
end

function homodyne_q!(Γ::State, m::Int, ::NoNoise)
    q = Γ.μ[2m-1]
    vacuum_project!(Γ, m)
    return q
end

function photon_number(Γ::State)
    return 1/2 * (norm(Γ.μ)^2 + sum(Γ.Σ[i]-1/2 for i in diagind(Γ.Σ)))
end
