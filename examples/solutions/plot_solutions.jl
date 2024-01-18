# Script to plot graphs of various nodes solutions in this folder

using JSON
import LaplacianOpt as LOpt

num_nodes = 25
instance = 2

filename = joinpath(@__DIR__, "$(num_nodes)_nodes/$(num_nodes)_$(instance).json")

open(filename) do f
    json_str = read(f, String)
    dict = JSON.parse(json_str)

    adjacency_matrix = zeros(num_nodes, num_nodes)
    for k in keys(dict["graph_adjacency"])
        e = dict["graph_adjacency"][k]
        adjacency_matrix[e[1], e[2]] = 1
        adjacency_matrix[e[2], e[1]] = 1
    end

    return LOpt.plot_graphviz(
        adjacency_matrix,
        display_edge_weights = false,
        edge_length = 0.5,
    )
end
