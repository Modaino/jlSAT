using DifferentialEquations
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
        clauses[i] = [(C[m,i], (collect∘zip)(js,C[m,j] for j ∈ js)) for (js,m) ∈ zip(j_inds,m_inds)]
    end
    return [clauses...]
end

function G_rule!(du, u, p, t)
    # constants
    β = 0.3
    η = 1.0
    ζ = p[2]
    c = p[1]
    M, N = size(c)

    # variables
    q = u[1:N]
    σ_sq = u[N+1:2N]
    e = u[2*N+1:3N]
    
    for i in 1:N
        for k in 1:N
            # γ_ik 
            du[i,k] = e[i] * sum(0.125*c[j,i]*c[j,k]*ζ* prod(l == i || l==k ? 1 : (1-ζ*c[j,l]*q[l]) for l in 1:N)/(2*sqrt(η)) for j in 1:M)
            # Γ_ik
            du[i,2N+k] = -β*e[i] * ζ^2 *q[i]/sqrt(η)
        end
        # γ_ik 
        du[i,i] += 2*sqrt(η)*(σ_sq[i]-0.5)
    end
end

function G_rule3!(dG, u, params, t)
    clauses, β, η, ς, p = params
    q, σ², e = u.x
    N = length(q)
    dG .= 0
    for i ∈ eachindex(q)
        dG[i,i] = 2*sqrt(η) * (σ²[i]-0.5)
        dG[i,2N+1:end] .= -β * e[i] * ς^2/sqrt(η) * q[i]
        for (Cmi,Cm) ∈ clauses[i]
            for (j,Cmj) ∈ Cm
                temp = 1
                for (k,Cmk) ∈ Cm
                    if k ≠ j && k ≠ i
                        temp *= 1 - ς * Cmk * q[k]
                    end
                end
                dG[i,j] += temp * (1/8) * ς / 2sqrt(η) * e[i] * Cmi * Cmj
            end
        end
        for (Cmi,Cm) ∈ clauses[i]
            temp = 1
            for (j,Cmj) ∈ Cm
                temp *= 1 - ς * Cmj * q[j]
            end
            dG[i,i] += (1/8) * ς / 2sqrt(η) * e[i] * temp
        end
    end
end

C = load_cnf(joinpath("SAT_problems","random3SATn50a3.42.cnf"))
M, N = size(C)
u = randn(3N)

# Base version
dG = zeros(N,3N)
@btime G_rule!(dG, u, (C,0.01), 0.)

# Sparse version
clauses = process_clauses(C)
params = (clauses, 0.3, 1.0, 0.01, 0.9)
u3 = ArrayPartition(u[1:N], u[N+1:2N], u[2N+1:3N])

dG3 = zeros(N,3N)
@btime G_rule3!(dG3, u3, params, 0.)
@test all(dG3 .≈ dG);
