using DifferentialEquations
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

function CTDS_rule!(du, u, c, t)
    M, N = size(c)
    s = u[1:N]
    a = u[N+1:end]
    for i in 1:N 
        du[i] = sum([2*a[m]*c[m,i]*(0.125*prod([ j!=i ? (1-c[m, j]*s[j]) : 1 for j in 1:length(s)]))^2*(1-c[m, i]*s[i]) for m in 1:M])
    end
    for m in 1:M
        du[N+m] = a[m]*0.125*prod([ (1-c[m, j]*s[j]) for j in 1:N])
    end
    return nothing
end

c = load_cnf("SAT_problems\\easy.cnf")
M, N = size(c)
u0 = ones(N+M)

u0[1:N] = [-0.048918786692759286, 0.49505812133072413, 0.8455355137115523, -0.373167063744003]

tspan = (0.0,100.0)
prob = ODEProblem(CTDS_rule!, u0, tspan, c)
sol = solve(prob)
plot(sol, xaxis="t", vars=(1:4), label=["s1" "s2" "s3" "s4" ])

