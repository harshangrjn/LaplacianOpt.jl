function visualize_solution(
    results::Dict{String,Any},
    data::Dict{String,Any};
    visualizing_tool = "tikz",
    plot_file_format = "pdf",
    display_edge_weights = true,
)
    adjacency_full_graph = data["adjacency_augment_graph"]
    (data["num_edges_existing"] > 0) &&
        (adjacency_full_graph += data["adjacency_base_graph"])

    if results["solution_type"] == "optimal"
        solution_mat = abs.(results["solution"]["z_var"])
    elseif results["solution_type"] == "heuristic"
        solution_mat = abs.(results["heuristic_solution"])
    end

    if sum(
        (isapprox.(solution_mat, 0, atol = 1E-5)) +
        (isapprox.(solution_mat, 1, atol = 1E-5)),
    ) == data["num_nodes"]^2
        Memento.info(_LOGGER, "Plotting the graph of integral solution")
        if visualizing_tool == "tikz"
            LOpt.plot_tikzgraph(
                solution_mat .* adjacency_full_graph,
                plot_file_format = plot_file_format,
                display_edge_weights = display_edge_weights,
            )

        elseif visualizing_tool == "graphviz"
            LOpt.plot_graphviz(
                solution_mat .* adjacency_full_graph,
                display_edge_weights = display_edge_weights,
            )
        end
    else
        Memento.info(
            _LOGGER,
            "Cannot plot as the obtained solutions are non-integral; fractional values can be found in the results dictionary",
        )
        return
    end
end

function plot_tikzgraph(
    adjacency_matrix::Matrix{<:Number};
    plot_file_format = "pdf",
    display_edge_weights = false,
)
    num_nodes = size(adjacency_matrix)[1]
    solution_graph = Graphs.SimpleGraph(num_nodes)
    edge_labels = Dict{Tuple{Int64,Int64},String}()

    for i in 1:(num_nodes-1)
        for j in (i+1):num_nodes
            if !isapprox(abs(adjacency_matrix[i, j]), 0, atol = 1E-6)
                Graphs.add_edge!(solution_graph, i, j)
                edge_labels[(i, j)] = string(ceil(adjacency_matrix[i, j], digits = 2))
            end
        end
    end

    if display_edge_weights
        # Plot with edge weights - graphs do not look great
        t = TikzGraphs.plot(
            solution_graph,
            edge_labels = edge_labels,
            node_style = "draw, rounded corners, fill=blue!10",
            TikzGraphs.Layouts.SpringElectrical(),
        )
    else
        t = TikzGraphs.plot(
            solution_graph,
            node_style = "draw, rounded corners, fill=blue!10",
            TikzGraphs.Layouts.SpringElectrical(),
        )
    end

    file_path =
        joinpath(dirname(pathof(LaplacianOpt)), "..", "examples/plots/plot_$(num_nodes)")

    if plot_file_format == "pdf"
        TikzPictures.save(TikzPictures.PDF(file_path), t)
    elseif plot_file_format == "tex"
        TikzPictures.save(TikzPictures.TEX(file_path), t)
    end
end

function plot_graphviz(
    adjacency_matrix::Matrix{<:Number};
    display_edge_weights = true,
    edge_length = 0.5,
)
    num_nodes = size(adjacency_matrix)[1]
    file_path = joinpath(
        dirname(pathof(LaplacianOpt)),
        "..",
        "examples/plots/plot_$(num_nodes).dot",
    )

    open(file_path, "w") do file
        write(file, "graph G { \n")
        write(file, "layout=neato; \n")
        write(file, "size=\"10,5\"; \n")
        # write(file, "rankdir=LR; \n")
        write(
            file,
            "node [fontname=\"Helvetica\", fontsize=7, shape = circle, width=0.15, fixedsize=true, style=\"filled\", fillcolor=\"cyan\"]; \n",
        )

        for i in 1:(num_nodes-1)
            for j in (i+1):num_nodes
                if !isapprox(abs(adjacency_matrix[i, j]), 0, atol = 1E-6)
                    if display_edge_weights
                        w_ij = string(ceil(adjacency_matrix[i, j], digits = 3))
                        write(
                            file,
                            "$i -- $j [label = \"$w_ij\", fontsize=9, fontname=\"Helvetica\", len=$edge_length]; \n",
                        )
                    else
                        write(
                            file,
                            "$i -- $j [fontsize=9, fontname=\"Helvetica\", len=$edge_length]; \n",
                        )
                    end
                end
            end
        end

        return write(file, "}")
    end

    # Further to convert the .dot in to a pdf file, use the following command in the terminal 
    # execute the dot_to_pdf.sh script file in the examples/plots folder

end
