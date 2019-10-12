function generate_graphs(verts = 4)
    graphs = Array{DiGraph}([])
    for _ in 1 : 10^4
        mat = rand(0 : 1, (verts, verts))
        for i in 1 : verts
            mat[i, i] = 0
        end
        gmat = DiGraph(mat)
        no_isom = true
        for g in graphs
            if has_isomorph(g, gmat)
                no_isom = false
            end 
        end
        if no_isom
            push!(graphs, gmat)
        end
    end 
    return graphs
end


function printed()
    out = open("GRAPH4.txt", "w")
    println(out, "[")
    for i in 1 : length(gs_sorted)
        mat = Matrix(adjacency_matrix(gs_sorted[i]))
        println(out, "[")
        for i in 1 : 4
            for j in 1 : 4
                print(out, mat[i, j])
                print(out, " ")
            end
            println(out, ";")
        end
        println(out, "],")
    end
    println(out, "]")
    close(out)
end


