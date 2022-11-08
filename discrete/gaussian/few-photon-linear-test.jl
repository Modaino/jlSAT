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
Δt = 0.01    # Continuous to discrete parameter
κ = 2.0      # Outcoupling rate
γ_cac = 5.0  # CAC gain rate
n_cac = 2.0  # CAC target amplitude-squared

# Compute discrete-time parameters
R_loss = 1-(1-Δt)^2
R_out = 2Δt * κ
G = γ_cac * Δt / 2n_cac
A = 2n_cac
# T_decay = -1/log(sqrt(1-R_out)*sqrt(1-R_loss))

# No proportional gain for now
# T_prop = 1.5*T_decay / min(40,max(1, 15*log10(A)+10.5)) * A

# Set simulation parameters
T_sim = 20 * T_decay
N_traj = 4
seed = rand(UInt32)

# Load problem
C = load_cnf(joinpath(pwd(),"SAT_problems","random3SATn50a3.42.cnf"))
M, N = size(C)
problem = process_clauses(C)

# Run simulation
rng = MersenneTwister(seed)
opo = OPO(NoGainModel(), R_loss, R_out)
csm = CSM(opo, problem)
q, σ², w, e = simulate(csm, G, A, T_sim, N_traj, rng);

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
