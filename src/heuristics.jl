# -- UNDER CONSTRUCTION -- #
# Heuristics to maximize algebraic connectivity of weighted graphs
using Graphs
using LinearAlgebra


function heuristic_kopt(data::Dict{String,Any}, kopt_parameter::Int, num_of_central_nodes_verifier::Int)
    num_nodes = data["num_nodes"]
    adjacency_augment_graph = data["adjacency_augment_graph"]
    adjacency_base_graph = data["adjacency_base_graph"]
    num_edges_existing = data["num_edges_existing"]
    num_edges_to_augment = data["num_edges_to_augment"]
    augment_budget = data["augment_budget"]
    is_base_graph_connected = data["is_base_graph_connected"]

    # Making list of all possible edges to augment
    edge_augment_list = collect(edges(SimpleGraph(adjacency_augment_graph)))

    # Making list of tuples of edges and edge weights
    edge_edgeweight_list = Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge,Float64}}(undef,length(edge_augment_list),)

    for (index, edge) in enumerate(edge_augment_list)
        edge_edgeweight_list[index] = (edge, adjacency_augment_graph[src(edge), dst(edge)])
    end

    # Sorting the list of tuples in descending order of edge weights
    sorted_edge_edgeweight_list = sort(edge_edgeweight_list, by = x -> x[2], rev = true)
    #@show sorted_edge_edgeweight_list

    if num_edges_existing == 0
        if augment_budget  == num_nodes - 1 #spanning tree

            # Priority central nodes based on sum of edge weights
            priority_central_nodes_list = priority_central_nodes(adjacency_augment_graph, num_nodes) 
            #@show priority_central_nodes_list

            # To keep track of adjacency matrix of highest algebraic connectivity and it's algebraic connectivity
            adjacency_graph_list = Vector{Tuple{Any,Any}}(undef, num_of_central_nodes_verifier)

            Threads.@threads for index in 1:num_of_central_nodes_verifier
                # Builds a spanning tree
                adjacency_graph_list[index] = build_span_tree(
                    num_nodes, 
                    adjacency_augment_graph, 
                    sorted_edge_edgeweight_list, 
                    priority_central_nodes_list[index], 
                )

                adjacency_graph_list[index] = refinement_span_tree(
                    adjacency_augment_graph,
                    sorted_edge_edgeweight_list,
                    adjacency_graph_list[index],
                    kopt_parameter,
                )
            end

            adjacency_graph_ac_tuple_star = adjacency_graph_list[sortperm(adjacency_graph_list, by = x -> x[2], rev = true)[1]]
            @show adjacency_graph_ac_tuple_star

        elseif augment_budget >= num_nodes 
            print("Under construction")
        end
    else
        print("Under construction")
    end

    return
end


function build_span_tree(
    num_nodes,
    adjacency_augment_graph,
    sorted_edge_edgeweight_list,
    central_node,
    )

    G = SimpleGraph(zeros(num_nodes, num_nodes)) # Starting graph
    uncon_nodes_set = Int64[] # Set contains all unconnected nodes
    uncon_nodes_set = collect(1:num_nodes)
    dir_con_set = Int64[] # Set contains nodes which are directly connected to central node
    sec_con_set = Int64[] # Set contains nodes which are connected to nodes in dir_con_set

    for j in eachindex(sorted_edge_edgeweight_list)
        if (src(sorted_edge_edgeweight_list[j][1]) == central_node || dst(sorted_edge_edgeweight_list[j][1]) == central_node)
            
            # First edge of the graph
            add_edge!(G, src(sorted_edge_edgeweight_list[j][1]), dst(sorted_edge_edgeweight_list[j][1])) 

            # Adding to directly connected set
            src(sorted_edge_edgeweight_list[j][1]) == central_node && push!(dir_con_set, dst(sorted_edge_edgeweight_list[j][1])) 
            dst(sorted_edge_edgeweight_list[j][1]) == central_node && push!(dir_con_set, src(sorted_edge_edgeweight_list[j][1])) 

            #  Removing from unconnected set
            deleteat!(uncon_nodes_set, uncon_nodes_set .== src(sorted_edge_edgeweight_list[j][1])) 
            deleteat!(uncon_nodes_set, uncon_nodes_set .== dst(sorted_edge_edgeweight_list[j][1]))
            break
        end
    end
    #@show connected_components(G)

    while !is_connected(G) # Connecting nodes until a spanning tree is formed
        algebraic_connectivity_list = Vector{Float64}(undef, size(dir_con_set)[1] + 1)
        algebraic_connectivity_ranker = Int64[]

        for j in 1:(size(dir_con_set)[1]+1)
            for k in eachindex(sorted_edge_edgeweight_list)
                if j == 1
                    if (src(sorted_edge_edgeweight_list[k][1]) == central_node) && issubset([dst(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set) ||
                    (dst(sorted_edge_edgeweight_list[k][1]) == central_node) && issubset([src(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set)
                        
                        # Adding edge to central node
                        add_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))
                        
                        # Storing current algebraic connectivity
                        algebraic_connectivity_list[j] = algebraic_connectivity((weighted_adj_matrix(G, adjacency_augment_graph, num_nodes - size(uncon_nodes_set)[1] + 1)))
                        

                        # Removing the added edge
                        rem_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))
                        break
                    end
                else
                    if (src(sorted_edge_edgeweight_list[k][1]) == dir_con_set[j-1]) && issubset([dst(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set) ||
                    (dst(sorted_edge_edgeweight_list[k][1]) == dir_con_set[j-1]) && issubset([src(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set)
                        
                        # Adding edge to directly connected node
                        add_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))

                        # Storing current algebraic connectivity
                        algebraic_connectivity_list[j] = algebraic_connectivity((weighted_adj_matrix(G, adjacency_augment_graph, num_nodes - size(uncon_nodes_set)[1] + 1)))
                        
                        # Removing the added edge
                        rem_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))
                        break
                    end
                end
            end
        end

        #Sorting algebraic connectivities to choose the best incumbent edge
        algebraic_connectivity_ranker = sortperm(algebraic_connectivity_list, rev = true)[1]

        for k in eachindex(sorted_edge_edgeweight_list)
            if algebraic_connectivity_ranker == 1
                if (src(sorted_edge_edgeweight_list[k][1]) == central_node) && issubset([dst(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set)

                    # Adding edge to central node
                    add_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))

                    # Adding node to directly connected set
                    push!(dir_con_set, dst(sorted_edge_edgeweight_list[k][1]))
                    
                    # Removing node from unconnected nodes set
                    deleteat!(uncon_nodes_set, uncon_nodes_set .== dst(sorted_edge_edgeweight_list[k][1]))
                    break

                elseif (dst(sorted_edge_edgeweight_list[k][1]) == central_node) && issubset([src(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set)

                    # Adding edge to central node
                    add_edge!(G, dst(sorted_edge_edgeweight_list[k][1]), src(sorted_edge_edgeweight_list[k][1]))
                    
                    # Adding node to directly connected set
                    push!(dir_con_set, src(sorted_edge_edgeweight_list[k][1]))
                    
                    # Removing node from unconnected nodes set
                    deleteat!(uncon_nodes_set, uncon_nodes_set .== src(sorted_edge_edgeweight_list[k][1]))
                    break
                end
            else
                if (src(sorted_edge_edgeweight_list[k][1]) == dir_con_set[algebraic_connectivity_ranker - 1]) && issubset([dst(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set)

                    # Adding edge to directly connected node
                    add_edge!(G, dir_con_set[algebraic_connectivity_ranker - 1], dst(sorted_edge_edgeweight_list[k][1]))

                    # Adding node to secondary connected set
                    push!(sec_con_set, dst(sorted_edge_edgeweight_list[k][1]))

                    # Removing node from unconnected nodes set
                    deleteat!(uncon_nodes_set, uncon_nodes_set .== dst(sorted_edge_edgeweight_list[k][1]))
                    break

                elseif (dst(sorted_edge_edgeweight_list[k][1]) == dir_con_set[algebraic_connectivity_ranker - 1]) && issubset([src(sorted_edge_edgeweight_list[k][1])], uncon_nodes_set)

                    # Adding edge to directly connected node
                    add_edge!(G, dir_con_set[algebraic_connectivity_ranker - 1], src(sorted_edge_edgeweight_list[k][1]))

                    # Adding node to secondary connected set
                    push!(sec_con_set, src(sorted_edge_edgeweight_list[k][1]))

                    # Removing node from unconnected nodes set
                    deleteat!(uncon_nodes_set, uncon_nodes_set .== src(sorted_edge_edgeweight_list[k][1]))
                    break
                end
            end
        end
        #@show connected_components(G)
    end
    #@show G
    return (Matrix(adjacency_matrix(G)), algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))))
end

function refinement_span_tree(
    adjacency_augment_graph,
    sorted_edge_edgeweight_list,
    adjacency_graph_ac_tuple,
    kopt_parameter,
    )

    G = SimpleGraph(adjacency_graph_ac_tuple[1])
    # Refinement of Algebraic connectivity
    if kopt_parameter == 1
        cycle_basis_list = Vector{Int64}
        for i in eachindex(sorted_edge_edgeweight_list)
            if add_edge!(G, src(sorted_edge_edgeweight_list[i][1]), dst(sorted_edge_edgeweight_list[i][1]))
                if is_cyclic(G)
                    cycle_basis_list = cycle_basis(G)[1]
                    algebraic_connectivity_list = Vector{Float64}(undef, 0)
                    vertices_list = Vector{Any}(undef, 0)
                    for j = 1:(length(cycle_basis_list)-1)
                        for k = (j+1):length(cycle_basis_list)
                            if rem_edge!(G, cycle_basis_list[j], cycle_basis_list[k])
                                push!(algebraic_connectivity_list, algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))))
                                push!(vertices_list, (cycle_basis_list[j], cycle_basis_list[k]))
                                add_edge!(G, cycle_basis_list[j], cycle_basis_list[k])
                            end
                        end
                    end
                    rem_edge!(G, vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][1], vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][2])
                end
            end
        end
        updated_adjacency_graph_ac_tuple = update(G, adjacency_augment_graph, adjacency_graph_ac_tuple)

    elseif kopt_parameter == 2
        combinations = edge_combinations(length(sorted_edge_edgeweight_list), kopt_parameter)
        cycle_basis_list = Vector{Int64}
        for i in eachindex(combinations)
            if adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1])] == 0 && 
            adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1])] == 0
                add_edge!(G, src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1]))
                add_edge!(G, src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1]))
                if is_cyclic(G)
                    cycle_basis_list = cycle_basis(G)
                    algebraic_connectivity_list = Vector{Float64}(undef, 0)
                    vertices_list = Vector{Any}(undef, 0)
                    for j in 1:(length(cycle_basis_list[1])-1)
                        for k in (j+1):length(cycle_basis_list[1])
                            for l in 1:(length(cycle_basis_list[2])-1)
                                for m in (l+1):length(cycle_basis_list[2])
                                    if (cycle_basis_list[1][j], cycle_basis_list[1][k]) !== (cycle_basis_list[2][l], cycle_basis_list[2][m])
                                        if adjacency_matrix(G)[cycle_basis_list[1][j], cycle_basis_list[1][k]] == 1 && adjacency_matrix(G)[cycle_basis_list[2][l], cycle_basis_list[2][m]] == 1
                                            rem_edge!(G, cycle_basis_list[1][j], cycle_basis_list[1][k])
                                            rem_edge!(G, cycle_basis_list[2][l], cycle_basis_list[2][m])
                                            push!(algebraic_connectivity_list, algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))))
                                            push!(vertices_list, ((cycle_basis_list[1][j], cycle_basis_list[1][k]), (cycle_basis_list[2][l], cycle_basis_list[2][m])))
                                            add_edge!(G, cycle_basis_list[1][j], cycle_basis_list[1][k])
                                            add_edge!(G, cycle_basis_list[2][l], cycle_basis_list[2][m])
                                        end
                                    end
                                end
                            end
                        end
                    end
                    rem_edge!(G, vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][1][1], vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][1][2])
                    rem_edge!(G, vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][2][1], vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][2][2])
                end
            end
        end
        updated_adjacency_graph_ac_tuple = update(G, adjacency_augment_graph, adjacency_graph_ac_tuple)
        
    elseif kopt_parameter == 3
        combinations = edge_combinations(length(sorted_edge_edgeweight_list), kopt_parameter)
        cycle_basis_list = Vector{Int64}
        for i in eachindex(combinations)
            if adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1])] == 0 &&  
            adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1])] == 0 &&  
            adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][3]][1]), dst(sorted_edge_edgeweight_list[combinations[i][3]][1])] == 0
                add_edge!(G, src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1]))
                add_edge!(G, src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1]))
                add_edge!(G, src(sorted_edge_edgeweight_list[combinations[i][3]][1]), dst(sorted_edge_edgeweight_list[combinations[i][3]][1]))
                if is_cyclic(G)
                    cycle_basis_list = cycle_basis(G)
                    algebraic_connectivity_list = Vector{Float64}(undef, 0)
                    vertices_list = Vector{Any}(undef, 0)
                    for j in 1:(length(cycle_basis_list[1])-1)
                        for k in (j+1):length(cycle_basis_list[1])
                            for l in 1:(length(cycle_basis_list[2])-1)
                                for m in (l+1):length(cycle_basis_list[2])
                                    for p in 1:(length(cycle_basis_list[3])-1)
                                        for q in (p+1):length(cycle_basis_list[3])
                                            if (cycle_basis_list[1][j], cycle_basis_list[1][k]) !== (cycle_basis_list[2][l], cycle_basis_list[2][m]) && 
                                            (cycle_basis_list[1][j], cycle_basis_list[1][k]) !== (cycle_basis_list[3][p], cycle_basis_list[3][q]) && 
                                            (cycle_basis_list[3][p], cycle_basis_list[3][q]) !== (cycle_basis_list[2][l], cycle_basis_list[2][m])
                                                if adjacency_matrix(G)[cycle_basis_list[1][j], cycle_basis_list[1][k]] == 1 && adjacency_matrix(G)[cycle_basis_list[2][l], cycle_basis_list[2][m]] == 1 && adjacency_matrix(G)[cycle_basis_list[3][p], cycle_basis_list[3][q]] == 1
                                                    rem_edge!(G, cycle_basis_list[1][j], cycle_basis_list[1][k])
                                                    rem_edge!(G, cycle_basis_list[2][l], cycle_basis_list[2][m])
                                                    rem_edge!(G, cycle_basis_list[3][p], cycle_basis_list[3][q])
                                                    push!(algebraic_connectivity_list, algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))))
                                                    push!(vertices_list, ((cycle_basis_list[1][j], cycle_basis_list[1][k]), (cycle_basis_list[2][l], cycle_basis_list[2][m]), (cycle_basis_list[3][p], cycle_basis_list[3][q])))
                                                    add_edge!(G, cycle_basis_list[1][j], cycle_basis_list[1][k])
                                                    add_edge!(G, cycle_basis_list[2][l], cycle_basis_list[2][m])
                                                    add_edge!(G, cycle_basis_list[3][p], cycle_basis_list[3][q])
                                                end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                    rem_edge!(G, vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][1][1], vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][1][2])
                    rem_edge!(G, vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][2][1], vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][2][2])
                    rem_edge!(G, vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][3][1], vertices_list[sortperm(algebraic_connectivity_list, rev=true)[1]][3][2])
                end
            end
        end
        updated_adjacency_graph_ac_tuple = update(G, adjacency_augment_graph, adjacency_graph_ac_tuple)
    end
    return updated_adjacency_graph_ac_tuple
end