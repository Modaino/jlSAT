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

function is_satisfied(spin_config, c)
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

function SATCFC_rule!(du, u, c, t)
    # parameters
    p = 0.001
    beta = 0.15
    M, N = size(c)
    
    # dummy variables
    x = u[1:N]
    e = u[N+1:2N]
    f = zeros(N)
    
    # updating derivatives
    for i in 1:N
        f[i] = e[i]*sum([c[j,i]*(0.125*prod([ k!=i ? (1-c[j, k]*x[k]) : 1 for k in 1:N])) for j in 1:M])
        du[i] = (p-1)*x[i]-x[i]^3+f[i]
        du[i+N] = -beta*e[i]*(f[i]^2-(p-2)^2)
    end
end

# Load SAT problem
c = load_cnf("SAT_problems\\random3SATn50a3.42.cnf")
M, N = size(c)

# initial condition(s)
u0 = 2*rand(2*N)-ones(2*N)
u0[N+1:2*N] = ones(N)

# terminate integration on solution
condition(u, t, integrator) = is_satisfied([s>0 ? 1 : -1 for s in u], c)
affect!(integrator) = terminate!(integrator)
cb = DiscreteCallback(condition, affect!)

# define problem and run simulation
tspan = (0.0,50.0)
prob = ODEProblem(SATCFC_rule!, u0, tspan, c)
sol = solve(prob, Tsit5(), callback = cb)

# Plotting and results
p1 = plot(sol, xaxis="t", yaxis="spins", vars=(1:N))
p2 = plot(sol, xaxis="t", yaxis="errors", vars=(N+1:2*N))
to_be_plotted = zeros(length(sol.t))
for idx in 1:length(sol.t)
    to_be_plotted[idx] = number_of_satisfied_clauses([sol[idx][i] > 0 ? 1 : -1 for i in 1:N], c)
end
p3 = plot(sol.t, 100*to_be_plotted/M, xaxis="analog time", yaxis="Satisfied clauses [%]")

plot(p1, p2, p3, layout= (3, 1), legend = false)

