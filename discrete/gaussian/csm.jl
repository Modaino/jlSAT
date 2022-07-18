import Random: AbstractRNG

struct CSM{T<:OPO,C}
    opo::T
    problem::C
    Γ::Vector{State}
    w::Vector{Float64}
    e::Vector{Float64}
    function CSM(opo::OPO, problem)
        N = length(problem)
        new{typeof(opo),typeof(problem)}(opo, problem, [Vacuum(2) for i ∈ 1:N], zeros(N), ones(N))
    end
end

function reset!(csm::CSM)
    for Γk ∈ csm.Γ
        vacuum_project!(Γk, 1); vacuum_project!(Γk, 2)
    end
    csm.w .= 0
    csm.e .= 1
end

function roundtrip!(csm::CSM, β::Real, G::Real, A::Real, rng::AbstractRNG=Random.GLOBAL_RNG; opo_kwargs...)
    csm.e .+= -G .* (csm.w.^2 .- A) .* csm.e
    a = sqrt(A)
    for i ∈ eachindex(csm.Γ)
        fi = 0
        for (Cmi,Cm) ∈ problem[i]
            fim = (1/4) * Cmi
            for (j,Cmj) ∈ Cm
                fim *= 1 - Cmj * csm.w[j] / a
            end
            fi += fim
        end
        fi *= csm.e[i]
        displace_q!(csm.Γ[i], 1, fi)
        csm.w[i] = circulate!(csm.opo, β, csm.Γ[i], rng; opo_kwargs...)
    end
    return csm
end
