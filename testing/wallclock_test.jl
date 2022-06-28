using DifferentialEquations
# using StaticArrays
using BenchmarkTools

function load_cnf(file_name)
    c = Nothing
    open(file_name) do file
        for (idx, line) in enumerate(eachline(file))
            if idx == 1
                N = parse(Int32, split(line, " ")[3])
                M = parse(Int32, split(line, " ")[4])
                c = zeros(M,N)
            else
                variables = split(line, " ")
                for var_str in variables
                    var = parse(Int32, var_str)
                    if var != 0
                        if var > 0
                            c[idx-1, var] = 1
                        elseif var < 0
                            c[idx-1, -var] = -1
                        end
                    end
                end
            end
        end
    end
    return c
end

function process_clauses(C)
    clauses = Vector{Any}(undef, N)
    for (i,Ci) ∈ enumerate(eachcol(C))
        m_inds = findall(Ci .≠ 0)
        j_inds = [filter!(j->j≠i, findall(C[m,:].≠0)) for m ∈ m_inds]
        clauses[i] = [(C[m,i], (collect∘zip)(js,C[m,j] for j ∈ js)) for (js,m) ∈ zip(j_inds,m_inds)]
        # clauses[i] = [(C[m,i], SVector((collect∘zip)(js,C[m,j] for j ∈ js)...)) for (js,m) ∈ zip(j_inds,m_inds)]
    end
    return [clauses...]
    # return SVector(clauses...)
end

function F_rule3!(du, u, params, t)
    clauses, β, η, ς, p = params
    q, σ², e = u.x
    dq, dσ², de = du.x
    for i ∈ eachindex(q)
        # The following line is nice but type-unstable for some reason
        # ςfi = e[i] * sum(Cni/8 * prod(1-ς*Cnj*q[j] for (j,Cnj) ∈ Cn) for (Cni,Cn) ∈ clauses[i])
        ςfi = 0.0
        for (Cmi,Cm) ∈ clauses[i]
            fim = (1/8) * Cmi
            for (j,Cmj) ∈ Cm
                fim *= 1-ς*Cmj*q[j]
            end
            ςfi += fim
        end
        ςfi *= e[i]
        dq[i] = (p-1) * q[i] - ς^2 * q[i]^3 + ςfi / ς
        dσ²[i] = 2p * σ²[i] - 2 * (σ²[i]-0.5) - 4η * (σ²[i]-0.5)^2 - (ς*q[i])^2 * (3σ²[i]-1)
        de[i] = -β * e[i] * (ςfi^2 - (p-2)^2)
    end
end

C = load_cnf(joinpath("SAT_problems","random3SATn15a4.266666666666667.cnf"))
M, N = size(C)
clauses = process_clauses(C)

params3 = (clauses, 0.3, 1.0, 0.01, 0.9)

u0 = 2*rand(3*N)-ones(3*N)
u0[N+1:3*N] = ones(2N)

u = ArrayPartition(u0[1:N], u0[N+1:2N], u0[2N+1:3N])
du = zero(u)

tspan = (0.0,50.0)
prob3 = ODEProblem(F_rule3!, u, tspan, params3)
# @btime solve(prob, Tsit5())
# @btime F_rule3!(du, u, params3, 0)


function F_rule!(du, u, p, t)
    # constants
    β = 0.3
    η = 1.0
    ζ = p[2]
    r = 0.9
    c = p[1]
    M, N = size(c)

    q = u[1:N]
    σ_sq = u[N+1:2N]
    e = u[2*N+1:3N]
    f = zeros(N)
    for i in 1:N
        f[i] = (1/ζ)* e[i]*sum(c[j,i]*(0.125*prod( k!=i ? (1-ζ*c[j, k]*q[k]) : 1 for k in 1:N)) for j in 1:M)
        du[i] = (r-1)*q[i]-(ζ^2)*q[i]^3+f[i]
        du[i+N] = 2*r*σ_sq[i] - 2*(σ_sq[i]-1/2) - 4*η*(σ_sq[i]-1/2)^2-2*(q[i]^2)*(ζ^2)*(3/2*σ_sq[i]-1/2)
        du[i+2N] = -β*e[i]*((ζ^2)*f[i]^2-(r-2)^2)
    end
end

params = (C, 0.01)
prob = ODEProblem(F_rule!, u, tspan, params)
# @btime solve(prob, Tsit5())
# @btime F_rule!(du, u, params, 0)
