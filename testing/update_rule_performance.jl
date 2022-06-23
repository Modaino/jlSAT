using DifferentialEquations
using StaticArrays
using BenchmarkTools
using Test

function lastLine(file)
    ea = eachline(file)
    local line
    for l in ea
        line = l
    end
    return line
end

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
        clauses[i] = [(C[m,i], SVector((collect∘zip)(js,C[m,j] for j ∈ js)...)) for (js,m) ∈ zip(j_inds,m_inds)]
    end
    return SVector(clauses...)
end

function F_rule!(du, u, p, t)
    c = p[1]
    # constants
    β = 0.3
    η = 1.0
    ζ = p[2]
    p = 0.9
    M, N = size(c)
    q = u[1:N]
    σ_sq = u[N+1:2N]
    e = u[2*N+1:3N]
    f = zeros(N)
    for i in 1:N
        f[i] = (1/ζ)* e[i]*sum(c[j,i]*(0.125*prod( k!=i ? (1-ζ*c[j, k]*q[k]) : 1 for k in 1:N)) for j in 1:M)
        du[i] = (p-1)*q[i]-(ζ^2)*q[i]^3+f[i]
        du[i+N] = 2*p*σ_sq[i] - 2*(σ_sq[i]-1/2) - 4*η*(σ_sq[i]-1/2)^2-2*(q[i]^2)*(ζ^2)*(3/2*σ_sq[i]-1/2)
        du[i+2N] = -β*e[i]*((ζ^2)*f[i]^2-(p-2)^2)
    end
end

function F_rule1!(du, u, c, t)
    # constants
    β = 0.3
    η = 1.0
    ζ = 0.01#1e-1
    p = 0.9
    M, N = size(c)
    q = u[1:N]
    σ_sq = u[N+1:2N]
    e = u[2*N+1:3N]
    f = zeros(N)
    for i in 1:N
        f[i] = (1/ζ)* e[i]*sum(c[j,i]*(0.125*prod( k!=i ? (1-ζ*c[j, k]*q[k]) : 1 for k in 1:N)) for j in 1:M)
        du[i] = (p-1)*q[i]-(ζ^2)*q[i]^3+f[i]
        du[i+N] = 2*p*σ_sq[i] - 2*(σ_sq[i]-1/2) - 4*η*(σ_sq[i]-1/2)^2-2*(q[i]^2)*(ζ^2)*(3/2*σ_sq[i]-1/2)
        du[i+2N] = -β*e[i]*((ζ^2)*f[i]^2-(p-2)^2)
    end
end

function F_rule2!(du, u, c, t)
    # constants
    β = 0.3
    η = 1.0
    ζ = 0.01#1e-1
    p = 0.9
    M, N = size(c)
    q = @view u[1:N]
    σ_sq = @view u[N+1:2N]
    e = @view u[2*N+1:3N]
    f = zeros(N)
    for i in 1:N
        f[i] = (1/ζ)* e[i]*sum(c[j,i]*(0.125*prod( k!=i ? (1-ζ*c[j, k]*q[k]) : 1 for k in 1:N)) for j in 1:M)
        du[i] = (p-1)*q[i]-(ζ^2)*q[i]^3+f[i]
        du[i+N] = 2*p*σ_sq[i] - 2*(σ_sq[i]-1/2) - 4*η*(σ_sq[i]-1/2)^2-2*(q[i]^2)*(ζ^2)*(3/2*σ_sq[i]-1/2)
        du[i+2N] = -β*e[i]*((ζ^2)*f[i]^2-(p-2)^2)
    end
end

function F_rule3!(du, u, params, t)
    clauses, β, η, ς, p = params
    q, σ², e = u.x
    dq, dσ², de = du.x
    for i ∈ eachindex(q)
        # The following line is nice but type-unstable for some reason
        # ςfi = e[i] * sum(Cni/8 * prod(1-ς*Cnj*q[j] for (j,Cnj) ∈ Cn) for (Cni,Cn) ∈ clauses[i])
        ςfi = 0
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

C = load_cnf(joinpath("SAT_problems","random3SATn50a3.42.cnf"))
M, N = size(C)
u = randn(3N)

# Base version
du = zero(u)
@btime F_rule!(du, u, (C,0.01), 0.) # display(@benchmark F_rule!(du, u, (C,0.01), 0.))

# Version 1
du1 = zero(u)
@btime F_rule1!(du1, u, C, 0.)
@test all(du1 .≈ du)

# Version 2
du2 = zero(u)
@btime F_rule2!(du2, u, C, 0.)
@test all(du2 .≈ du)

# Version 3
clauses = process_clauses(C)
params = (clauses, 0.3, 1.0, 0.01, 0.9)
u3 = ArrayPartition(u[1:N], u[N+1:2N], u[2N+1:3N])

du3 = zero(u3)
@btime F_rule3!(du3, u3, params, 0.)
@test all(du3 .≈ du);
