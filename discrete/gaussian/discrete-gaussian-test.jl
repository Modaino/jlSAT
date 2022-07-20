using Polynomials
using Random
using PyPlot

include("noise.jl")
include("gaussian.jl")
include("crystal.jl")
include("opo.jl")
include("csm.jl")
include("sat.jl")

# Simulation code
function simulate(csm, β, G, A, T_sim, N_traj, rng; kwargs...)
    q = Vector{Vector{Vector{Float64}}}(undef, N_traj)
    σ² = deepcopy(q); w = deepcopy(q); e = deepcopy(q)
    for l ∈ 1 : N_traj
        reset!(csm)
        q[l] = Vector{Vector{Float64}}(undef, T_sim)
        σ²[l] = deepcopy(q[l]); w[l] = deepcopy(q[l]); e[l] = deepcopy(q[l])
        for t ∈ 1 : T_sim
            println("Traj $l, t = $t")
            roundtrip!(csm, β, G, A, rng; kwargs...)
            q[l][t] = [Γi.μ[1] for Γi in csm.Γ]
            σ²[l][t] = [Γi.Σ[1,1] for Γi in csm.Γ]
            w[l][t] = copy(csm.w)
            e[l][t] = copy(csm.e)
        end
    end
    return q, σ², w, e
end

# Set system parameters
T_decay = 16 # Cavity decay time, in roundtrips
η = 1         # Escape efficiency
r = 1.5       # Pump parameter
n_sat = 1e-1  # Saturation photon number
β_cac = 0.05   # CAC gain rate, in 1/T_decay
α_cac = 1.0   # CAC target photon number, in units of n_sat

# Set simulation parameters
T_sim = 20 * T_decay
N_traj = 1
seed = rand(UInt32)

# Calculate OPO parameters
R_out = η * (1 - exp(-2/T_decay))
R_loss = (1/η - 1) * R_out / (1 - R_out)
ϵτ_nl = sqrt(4 / T_decay / n_sat)
β = sqrt(2) * r / ϵτ_nl / T_decay

# Generate feedback parameters
n_cac = α_cac * n_sat
G = (β_cac/T_decay) / R_out / 2n_cac
κ = R_out/2; g = ϵτ_nl^2/4; p = r/T_decay
cfs = [p - n_cac*g/9, 2(p-(1+κ)-n_cac*g/9), -4κ]
δ = 0.5 + R_out * (roots∘Polynomial)(cfs)[2]
A = 2n_cac * R_out + δ

# Load problem
C = load_cnf(joinpath(pwd(),"SAT_problems","random3SATn50a3.42.cnf"))
M, N = size(C)
problem = process_clauses(C)

# Run simulation
rng = MersenneTwister(seed)
opo = OPO(GaussianGainModel(ϵτ_nl), R_loss, R_out)
csm = CSM(opo, problem)
q, σ², w, e = simulate(csm, β, G, A, T_sim, N_traj, rng);

# Plot results
fig, axs = plt.subplots(2, 2, figsize=(10,6))
axs[1,1].plot(q[1], alpha=0.3)
axs[1,2].plot(σ²[1], alpha=0.3)
axs[2,1].plot(e[1], alpha=0.3)
axs[2,2].plot(w[1], alpha=0.3)

axs[1,1].set_ylim(maximum(abs.(axs[1,1].get_ylim())) * [-1,1])
axs[1,2].set_ylim(0, axs[1,2].get_ylim()[2])
axs[2,1].set_ylim(0, axs[2,1].get_ylim()[2])
axs[2,2].set_ylim(maximum(abs.(axs[2,2].get_ylim())) * [-1,1])

axs[1,1].axhline(sqrt(α_cac * 2n_sat), linestyle="--", color=:k)
axs[1,1].axhline(-sqrt(α_cac * 2n_sat), linestyle="--", color=:k)

for i ∈ 1 : size(axs,1)
    for j ∈ 1 : size(axs,2)
        axs[i,j].set_xlim(0, T_sim)
        axs[i,j].tick_params(labelbottom=i==size(axs,1))
    end
end
