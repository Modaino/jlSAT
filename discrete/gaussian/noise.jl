import Random: AbstractRNG, randn
import Base: RefValue

struct NoNoise <: AbstractRNG end

struct ThermalNoise <: AbstractRNG
    σ_th::Float64
    rng::MersenneTwister
end

function randn(rng::ThermalNoise)
    return randn(rng.rng) * rng.σ_th
end

# Normally-distributed random downsampler, used to resample Wiener processes
struct NDRandnDownSampler <: AbstractRNG
    values::Vector{Float64}
    index::RefValue{Int64}
    function NDRandnDownSampler(noise::Matrix{Float64}, M::Int)
        T = div(size(noise,2), M)
        values = zeros(size(noise,1), T)
        for j in 0 : (T-1)
            for k in 1 : M
                values[:,j+1] .+= noise[:,j*M+k]
            end
        end
        values ./= sqrt(M)
        return new(vcat(values...), Ref(1))
    end
end

function randn(rng::NDRandnDownSampler)
    val = rng.values[rng.index[]]
    rng.index[] += 1
    return val
end
