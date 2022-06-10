using DifferentialEquations
using DynamicalSystems
using Plots

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

function number_of_satisfied_clauses(state, c)
    M, N = size(c)
    result = 0
    for m in 1:M
        satisfied = 0
        for j in 1:N
            if c[m,j] == state[j]
                satisfied = 1
                break
            end
        end
        result += satisfied
    end
    return result
end

function satisfied(spin_config, c)
    function check_clause(row, state)
        for (index,elem) in enumerate(row)
            if elem == state[index]
                return true
            end
        end
    end

    incorrect_flag = false
    for clause in eachrow(c)
        if check_clause(clause, spin_config) != true
            incorrect_flag = true
            break
        end
    end

    if incorrect_flag
        return false
    end
    return true
end

function G_rule!(du, u, c, t)
    # constants
    r = 1
    η = 1.0
    ζ = 1e-2
    β = 0.15# 2.5e-9
    M, N = size(c)

    # variables
    q = u[1:N]
    σ_sq = u[N+1:2N]
    e = u[2*N+1:3N]

    
    for i in 1:N
        # γ_ik 
        du[i,i] = 2*sqrt(η)*(σ_sq[i]-0.5)
        for k in 1:N
            # γ_ik 
            du[i,k] = e[i] * sum([ 0.125*c[j,i]*c[j,k]*ζ* prod([l == i || l==k ? 1 : (1-ζ*c[j,l]*q[l]) for l in 1:N])/(2*sqrt(η)) for j in 1:M])
            # Γ_ik
            du[i,2N+k] = -β*e[i] * ζ^2 *q[i]/sqrt(η)
        end
    end
end

function F_rule!(du, u, c, t)
    # constants
    r = 0.1
    η = 1.0
    ζ = 1e-2
    β = 0.15# 2.5e-9
    M, N = size(c)

    # variables
    q = u[1:N]
    σ_sq = u[N+1:2N]
    e = u[2*N+1:3N]

    # placeholders 
    f_noiseless = zeros(N)
    for i in 1:N
        # placeholders
        f_noiseless[i] = e[i]*sum([0.125*c[j,i]*prod([k==i ? 1 : 1-ζ*c[j,k]*q[k] for k in 1:N ]) for j in 1:M]) #squared term missing: only use for SAT withouth trivial clauses
        # updating rule
        du[i] = ((r-1)-(ζ^2)*q[i]^2)*q[i] + f_noiseless[i]
        du[i+N] = 2*r*σ_sq[i] - 2*(σ_sq[i]-1/2) - 4*η*(σ_sq[i]-1/2)^2-2*(ζ^2)*(q[i]^2)*(3/2*σ_sq[i]-1/2)
        du[i+2N]= -β*e[i]*ζ^2*(q[i]^2 - (r-2)^2 + 1/(4*η))
    end
end

c = load_cnf("SAT_problems\\random3SATn15a4.266666666666667.cnf")
M, N = size(c)

# initial & boundary conditions
tspan = (0.0, 500.0)
u0 = 0.5*ones(3N)
u0[1:N] = zeros(N) # spin variables

# define problem and run simulation
prob = SDEProblem(F_rule!, G_rule!, u0, tspan, c, noise_rate_prototype=zeros(3N, 3N))
sol = solve(prob,LambaEulerHeun())

# Plotting and results
p1 = plot(sol, xaxis="analog time", yaxis="q", vars=(1:N))
p2 = plot(sol, xaxis="analog time", yaxis="σ^2", vars=(N+1:2*N))
to_be_plotted = zeros(length(sol.t))
for idx in 1:length(sol.t)
    to_be_plotted[idx] =1/M * number_of_satisfied_clauses([sol[idx][i] > 0 ? 1 : -1 for i in 1:N], c)
end
p3 = plot(sol.t, to_be_plotted, xaxis="analog time", yaxis="Ratio of satisfied clauses")
p4 = plot(sol, xaxis="analog time", yaxis="errors", vars=(2*N+1:3*N))
plot(p1, p2, p3, p4, layout= (2, 2), legend = false)
