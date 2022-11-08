import Random: AbstractRNG

struct OPO{T<:GainModel}
    𝒳_xtal::T
    ℬ_loss::Symplectic
    ℬ_out::Symplectic
    function OPO(𝒳_xtal::GainModel, R_loss::Real, R_out::Real)
        ℬ_loss = Beamsplitter(1-sqrt(1-R_loss))
        ℬ_out = Beamsplitter(R_out)
        return new{typeof(𝒳_xtal)}(𝒳_xtal, ℬ_loss, ℬ_out)
    end
end

function circulate!(opo::OPO, β::Real, Γ::State, rng::AbstractRNG=Random.GLOBAL_RNG; λ_in::Real=1, de_kwargs...)
    vacuum_project!(Γ, 2); update!(opo.ℬ_loss, Γ)
    vacuum_project!(Γ, 2); propagate!(opo.𝒳_xtal, β, Γ; de_kwargs...)
    vacuum_project!(Γ, 2); update!(opo.ℬ_loss, Γ)
    squeeze_project!(Γ, 2, λ_in); update!(opo.ℬ_out, Γ)
    return homodyne_q!(Γ, 2, rng)
end

# Simplified for no crystal, no losses
function circulate!(opo::OPO{NoGainModel}, Γ::State, rng::AbstractRNG=Random.GLOBAL_RNG; λ_in::Real=1, de_kwargs...)
    squeeze_project!(Γ, 2, λ_in); update!(opo.ℬ_out, Γ)
    return homodyne_q!(Γ, 2, rng)
end

# More simplified for no squeezing at measurement
function circulate!(opo::OPO{NoGainModel}, Γ::State, rng::AbstractRNG=Random.GLOBAL_RNG; de_kwargs...)
    vacuum_project!(Γ, 2); update!(opo.ℬ_out, Γ)
    return homodyne_q!(Γ, 2, rng)
end
