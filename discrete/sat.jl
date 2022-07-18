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
