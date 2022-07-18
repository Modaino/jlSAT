import DifferentialEquations: ODEProblem, solve, ODEFunction

abstract type GainModel end

struct GaussianGainModel{T<:ODEProblem} <: GainModel
    ode::T
    function GaussianGainModel(ϵτ_nl::Real)
        eoms = ODEFunction(
            function(X′::Vector{<:Real}, X::Vector{<:Real}, p::Any, t::Real)
                q, u, δqδq, δpδp, δuδu, δvδv, δqδu, δpδv = X
                X′[1] = u * q + δqδu + δpδv
                X′[2] = -1/2 * (q^2 + δqδq - δpδp)
                X′[3] = 2 * (u * δqδq + q * δqδu)
                X′[4] = -2 * (u * δpδp - q * δpδv)
                X′[5] = -2 * q * δqδu
                X′[6] = -2 * q * δpδv
                X′[7] = u * δqδu + q * (δuδu - δqδq)
                X′[8] = -u * δpδv + q * (δvδv - δpδp)
            end
        )
        ode = ODEProblem(eoms, zeros(8), (0.0, ϵτ_nl/sqrt(2)))
        return new{typeof(ode)}(ode)
    end
end
function propagate!(𝒳::GaussianGainModel, β::Real, Γ::State; de_kwargs...)
    vacuum_project!(Γ, 2); displace_q!(Γ, 2, β)
    μ = Γ.μ; Σ = Γ.Σ; u0 = 𝒳.ode.u0
    u0[1], u0[2] = μ[1], μ[3]
    u0[3], u0[4], u0[5], u0[6] = Σ[1,1], Σ[2,2], Σ[3,3], Σ[4,4]
    u0[7], u0[8] = Σ[1,3], Σ[2,4]
    sol = solve(𝒳.ode; dense=false, save_everystep=false, de_kwargs...)
    μ[1], μ[3], Σ[1,1], Σ[2,2], Σ[3,3], Σ[4,4], Σ[1,3], Σ[2,4] = sol.u[end]
    Σ[3,1] = Σ[1,3]; Σ[4,2] = Σ[2,4]
    return Γ, sol
end

struct ClassicalGainModel{T<:ODEProblem} <: GainModel
    ode::T
    function ClassicalGainModel(ϵτ_nl::Real)
        eoms = ODEFunction(
            function(X′::Vector{<:Real}, X::Vector{<:Real}, p::Any, t::Real)
                X′[1] = X[2] * X[1]
                X′[2] = -0.5 * X[1]^2
            end
        )
        ode = ODEProblem(eoms, zeros(2), (0.0, ϵτ_nl/sqrt(2)))
        return new{typeof(ode)}(ode)
    end
end
function propagate!(𝒳::ClassicalGainModel, β::Real, Γ::State; de_kwargs...)
    vacuum_project!(Γ, 2); displace_q!(Γ, 2, β)
    𝒳.ode.u0[1], 𝒳.ode.u0[2] = Γ.μ[1], Γ.μ[3]
    sol = solve(𝒳.ode; dense=false, save_everystep=false, de_kwargs...)
    Γ.μ[1], Γ.μ[3] = sol.u[end]
    return Γ, sol
end

struct LinearGainModel <: GainModel end
function propagate!(::LinearGainModel, η::Real, Γ::State)
    Γ.μ[1] *= exp(η)
    Γ.Σ[1,1] *= exp(2*η)
    Γ.Σ[2,2] *= exp(-2*η)
end

struct NoGainModel <: GainModel end
propagate!(::NoGainModel, ::Real, Γ::State) = Γ

# Generic propagate! for saving intermediate points, dispatches to others
function propagate!(𝒳::GainModel, β::Real, Γ::State, Nsave::Int; de_kwargs...)
    Γ, sol = propagate!(𝒳, β, Γ; saveat=𝒳.ode.tspan[end]/Nsave, de_kwargs...)
    return sol.t, [[sol(t)[i] for t in sol.t] for i in eachindex(𝒳.ode.u0)]
end
