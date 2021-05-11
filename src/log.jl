function visualize_LOModel_solution(results::Dict{String, Any}, data::Dict{String, Any}; visualizing_tool = "tikz", plot_file_format = "pdf")
    
    if results["primal_status"] != MOI.FEASIBLE_POINT
        Memento.error(_LOGGER, "Non-feasible primal status. Graph solution may not be exact")
    end

    if !data["relax_integrality"]
        Memento.info(_LOGGER, "Plotting the graph of integral solution")
        
        if visualizing_tool == "tikz"
            LO.plot_tikzgraph(results["solution"]["z_var"], data["num_nodes"], data["instance"], data["edge_weights"], data["tol_zero"], plot_file_format)
        elseif visualizing_tool == "graphviz"
            LO.plot_graphviz(results["solution"]["z_var"], data["num_nodes"], data["instance"], data["edge_weights"], data["tol_zero"], plot_file_format)
        end

    else 
        Memento.info(_LOGGER, "Cannot plot as the obtained solutions are non-integral; fractional values can be found in the results dictionary")
        return
    end
    
end

function plot_tikzgraph(z_var::Matrix{Float64}, num_nodes::Int64, instance::Int64, edge_weights::Matrix{Float64}, tol_zero::Float64, plot_file_format::String)

    solution_graph = LG.SimpleGraph(num_nodes)

    edge_labels = Dict{Tuple{Int64, Int64}, String}()
    
    for i=1:(num_nodes-1)
        for j=(i+1):num_nodes

            if (1 - tol_zero) <= z_var[i,j] <= (1 + tol_zero)
                LG.add_edge!(solution_graph, i, j)
                edge_labels[(i,j)] = string(ceil(edge_weights[i,j], digits=2))
            end

        end
    end

    t = TikzGraphs.plot(solution_graph, node_style="draw, rounded corners, fill=blue!10", TikzGraphs.Layouts.SpringElectrical())

    # Plot with edge weights - graphs do not look great
    # t = TikzGraphs.plot(solution_graph, edge_labels=edge_labels, node_style="draw, rounded corners, fill=blue!10", TikzGraphs.Layouts.SpringElectrical())
    
    file_path = joinpath(dirname(pathof(LaplacianOpt)),"..", "examples/plots/plot_$(num_nodes)_$(instance)")

    if plot_file_format == "pdf"
        TikzPictures.save(TikzPictures.PDF(file_path), t)
    elseif plot_file_format == "tex"
        TikzPictures.save(TikzPictures.TEX(file_path), t)
    end
    
end

function plot_graphviz(z_var::Matrix{Float64}, num_nodes::Int64, instance::Int64, edge_weights::Matrix{Float64}, tol_zero::Float64, plot_file_format::String)
    
    file_path = joinpath(dirname(pathof(LaplacianOpt)),"..", "examples/plots/plot_$(num_nodes)_$(instance).dot")

    open(file_path, "w") do file

        write(file, "graph G { \n")
        write(file, "layout=neato; \n")
        write(file, "size=\"10,5\"; \n")
        # write(file, "rankdir=LR; \n")
        write(file, "node [fontname=\"Helvetica\", fontsize=20, shape = circle, width=0.4, fixedsize=true, style=\"filled\", fillcolor=\"0.650 0.200 1.000\"]; \n")

        for i=1:(num_nodes-1)
            for j=(i+1):num_nodes
    
                if (1 - tol_zero) <= z_var[i,j] <= (1 + tol_zero)
                    w_ij = string(ceil(edge_weights[i,j], digits=3))
                    write(file, "$i -- $j [label = \"$w_ij\", fontsize=9, fontname=\"Helvetica\"]; \n")
                end
    
            end
        end

        write(file, "}")
    end

    # Further to convert the .dot in to a pdf file, use the following command in the terminal 
    # execute the dot_to_pdf.sh script file in the examples/plots folder

end