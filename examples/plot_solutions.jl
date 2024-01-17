using JSON
import LaplacianOpt as LOpt

n = 25
instance = 2

filename = "solutions/$(n)_nodes/$(n)_$(instance).json"

open(filename) do f
    json_str = read(f, String)
    dict = JSON.parse(json_str)

    adjacency_matrix = zeros(n, n)
    for k in keys(dict["graph_adjacency"])
        e = dict["graph_adjacency"][k]
        @show e
        adjacency_matrix[e[1], e[2]] = 1
        adjacency_matrix[e[2], e[1]] = 1
    end

    return LOpt.plot_graphviz(adjacency_matrix, display_edge_weights = false)
end
