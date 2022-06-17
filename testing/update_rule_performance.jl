using BenchmarkTools
using Plots

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

# c = load_cnf("SAT_problems\\random3SATn60a4.267.cnf")

# M, N = size(c)
# ζ = 0.01

# @benchmark F_rule1!(zeros(3N), rand(3N), c, 16.1)
