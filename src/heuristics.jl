# Heuristics to maximize algebraic connectivity of weighted graphs
using Graphs
using LinearAlgebra


function heuristic_kopt(data::Dict{String,Any}, kopt_parameter::Int, num_central_nodes_verifier::Int, num_kopt_swaps_upperbound::Int)
    num_nodes = data["num_nodes"]
    adjacency_augment_graph = data["adjacency_augment_graph"]
    adjacency_base_graph = data["adjacency_base_graph"]
    num_edges_existing = data["num_edges_existing"]
    augment_budget = data["augment_budget"]
    is_base_graph_connected = data["is_base_graph_connected"]

    # Initializing the heuristic solutions
    adjacency_graph_star = Matrix{Float64}(zeros(num_nodes, num_nodes))
    algebraic_connectivity_star = 0

    # Making list of all possible edges to augment
    edge_augment_list = collect(edges(SimpleGraph(adjacency_augment_graph)))

    # Making list of tuples of edges and edge weights
    edge_edgeweight_list = Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge,Float64}}(undef,length(edge_augment_list))
    sorted_edge_edgeweight_list = Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge,Float64}}(undef,length(edge_augment_list))

    for (index, edge) in enumerate(edge_augment_list)
        edge_edgeweight_list[index] = (edge, adjacency_augment_graph[src(edge), dst(edge)])
    end

    # Sorting the list of tuples in descending order of edge weights
    sorted_edge_edgeweight_list = sort(edge_edgeweight_list, by = x -> x[2], rev = true)

    if num_edges_existing == 0
        if augment_budget  == num_nodes - 1 #spanning tree

            # Priority central nodes based on sum of edge weights
            priority_central_nodes_list = priority_central_nodes(adjacency_augment_graph, num_nodes) 

            # To keep track of adjacency matrix of highest algebraic connectivity and it's algebraic connectivity
            adjacency_graph_list = Vector{Tuple{Any,Any}}(undef, num_central_nodes_verifier)

            #JULIA_NUM_THREADS=auto #uncomment this line for multi threading (Requires atleast Julia 1.7)
            #Threads.@threads for index in 1:num_central_nodes_verifier #uncomment this line for multi threading
            for index in 1:num_central_nodes_verifier #comment this line for multi threading
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
                    num_kopt_swaps_upperbound,
                )
            end

            adjacency_graph_star, algebraic_connectivity_star = adjacency_graph_list[sortperm(adjacency_graph_list, by = x -> x[2], rev = true)[1]]

        elseif augment_budget >= num_nodes 
            print("Under construction")
        end
    else
        if is_base_graph_connected

            # Fiedler vector of base graph
            fiedler_vector_base_graph = fiedler_vector(adjacency_base_graph)

            # Making list of tuples of edges and fiedler weights (w_ij * (v_i - v_j)^2)
            edge_fiedlerweight_list = Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge,Float64}}(undef,length(edge_augment_list))
            sorted_edge_fiedlerweight_list = Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge,Float64}}(undef,length(edge_augment_list))

            for (index, edge) in enumerate(edge_augment_list)
                edge_fiedlerweight_list[index] = (edge, adjacency_augment_graph[src(edge), dst(edge)] * (fiedler_vector_base_graph[src(edge)] - fiedler_vector_base_graph[dst(edge)])^2)
            end

            # Sorting the list of tuples in descending order of fiedler weights
            sorted_edge_fiedlerweight_list = sort(edge_fiedlerweight_list, by = x -> x[2], rev = true)

            # Base graph
            G = SimpleGraph(adjacency_base_graph)

            # Adding edges based on fiedler weights to construct initial graph with required edges
            for i = 1:augment_budget
                add_edge!(G, src(sorted_edge_fiedlerweight_list[i][1]), dst(sorted_edge_fiedlerweight_list[i][1]))
            end

            adjacency_graph_star, algebraic_connectivity_star = refinement_tree(
                G,
                adjacency_base_graph,
                adjacency_augment_graph,
                sorted_edge_fiedlerweight_list,
                kopt_parameter,
                num_kopt_swaps_upperbound,
            )

        elseif !(is_base_graph_connected)
            print("Under construction")
        end
    end
    #  Implement the bridge to results here
    @show adjacency_graph_star, algebraic_connectivity_star
    return
end

function build_span_tree(
    num_nodes,
    adjacency_augment_graph,
    sorted_edge_edgeweight_list,
    central_node,
    )

    G = SimpleGraph(zeros(num_nodes, num_nodes)) # Starting graph
    uncon_set = Int64[] # Set contains all unconnected nodes
    uncon_set = collect(1:num_nodes)
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
            deleteat!(uncon_set, uncon_set .== src(sorted_edge_edgeweight_list[j][1])) 
            deleteat!(uncon_set, uncon_set .== dst(sorted_edge_edgeweight_list[j][1]))
            break
        end
    end

    while !is_connected(G) # Connecting nodes until a spanning tree is formed
        algebraic_connectivity_tracker = Vector{Tuple{Graphs.SimpleGraphs.SimpleEdge,Float64}}(undef,1)

        for j in 1:(size(dir_con_set)[1]+1)
            for k in eachindex(sorted_edge_edgeweight_list)
                if j == 1
                    if (src(sorted_edge_edgeweight_list[k][1]) == central_node) && issubset([dst(sorted_edge_edgeweight_list[k][1])], uncon_set) ||
                    (dst(sorted_edge_edgeweight_list[k][1]) == central_node) && issubset([src(sorted_edge_edgeweight_list[k][1])], uncon_set)
                        
                        # Adding edge to central node
                        add_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))
                        
                        # Storing current algebraic connectivity
                        algebraic_connectivity_tracker = ((src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1])),algebraic_connectivity((weighted_adjacency_matrix(G, adjacency_augment_graph, num_nodes - size(uncon_set)[1] + 1))))
                        
                        # Removing the added edge
                        rem_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))
                        break
                    end
                else
                    if (src(sorted_edge_edgeweight_list[k][1]) == dir_con_set[j-1]) && issubset([dst(sorted_edge_edgeweight_list[k][1])], uncon_set) ||
                    (dst(sorted_edge_edgeweight_list[k][1]) == dir_con_set[j-1]) && issubset([src(sorted_edge_edgeweight_list[k][1])], uncon_set)
                        
                        # Adding edge to directly connected node
                        add_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))

                        # Updating current algebraic connectivity
                        if algebraic_connectivity((weighted_adjacency_matrix(G, adjacency_augment_graph, num_nodes - size(uncon_set)[1] + 1))) > algebraic_connectivity_tracker[2]
                            algebraic_connectivity_tracker = ((src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1])),algebraic_connectivity((weighted_adjacency_matrix(G, adjacency_augment_graph, num_nodes - size(uncon_set)[1] + 1))))
                        end

                        # Removing the added edge
                        rem_edge!(G, src(sorted_edge_edgeweight_list[k][1]), dst(sorted_edge_edgeweight_list[k][1]))
                        break
                    end
                end
            end
        end

        # Adding edge
        add_edge!(G, algebraic_connectivity_tracker[1][1], algebraic_connectivity_tracker[1][2])
        
        if (algebraic_connectivity_tracker[1][1] ==  central_node) && issubset(algebraic_connectivity_tracker[1][2], uncon_set)
            # Adding node to directly connected set
            push!(dir_con_set, algebraic_connectivity_tracker[1][2])

            # Removing node from unconnected nodes set
            deleteat!(uncon_set, uncon_set .== algebraic_connectivity_tracker[1][2]) 

        elseif (algebraic_connectivity_tracker[1][2] ==  central_node) && issubset(algebraic_connectivity_tracker[1][1], uncon_set)
            # Adding node to directly connected set
            push!(dir_con_set, algebraic_connectivity_tracker[1][1])

            # Removing node from unconnected nodes set
            deleteat!(uncon_set, uncon_set .== algebraic_connectivity_tracker[1][1]) 

        elseif issubset(algebraic_connectivity_tracker[1][1], dir_con_set) && issubset(algebraic_connectivity_tracker[1][2], uncon_set)

            # Adding node to secondary connected set
            push!(sec_con_set, algebraic_connectivity_tracker[1][2])

            # Removing node from unconnected nodes set
            deleteat!(uncon_set, uncon_set .== algebraic_connectivity_tracker[1][2])

        elseif issubset(algebraic_connectivity_tracker[1][2], dir_con_set) && issubset(algebraic_connectivity_tracker[1][1], uncon_set)
            # Adding node to secondary connected set
            push!(sec_con_set, algebraic_connectivity_tracker[1][1])

            # Removing node from unconnected nodes set
            deleteat!(uncon_set, uncon_set .== algebraic_connectivity_tracker[1][1])
            
        end
        
    end
    return (Matrix(adjacency_matrix(G)), algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))))
end

function refinement_span_tree(
    adjacency_augment_graph,
    sorted_edge_edgeweight_list,
    adjacency_graph_ac_tuple,
    kopt_parameter,
    num_kopt_swaps_upperbound,
    )

    G = SimpleGraph(adjacency_graph_ac_tuple[1])
    # Refinement of Algebraic connectivity
    if kopt_parameter == 1
        cycle_basis_current = Vector{Int64}
        for i in eachindex(sorted_edge_edgeweight_list)
            if add_edge!(G, src(sorted_edge_edgeweight_list[i][1]), dst(sorted_edge_edgeweight_list[i][1]))
                if is_cyclic(G)
                    cycle_basis_current = cycle_basis(G)[1]
                    algebraic_connectivity_tracker = 0.0
                    vertices_tracker = [(undef, undef)]
                    for j = 1:(length(cycle_basis_current)-1)
                        for k = (j+1):length(cycle_basis_current)
                            if rem_edge!(G, cycle_basis_current[j], cycle_basis_current[k])
                                if algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))) > algebraic_connectivity_tracker
                                    algebraic_connectivity_tracker = algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G)))
                                    vertices_tracker = [(cycle_basis_current[j], cycle_basis_current[k])]
                                end
                                add_edge!(G, cycle_basis_current[j], cycle_basis_current[k])
                            end
                        end
                    end
                    rem_edge!(G, vertices_tracker[1][1], vertices_tracker[1][2])
                end
            end
        end
        updated_adjacency_graph_ac_tuple = update_kopt_adjacency!(G, adjacency_augment_graph, adjacency_graph_ac_tuple)

    elseif kopt_parameter == 2
        combinations = edge_combinations(length(sorted_edge_edgeweight_list), kopt_parameter)
        cycle_basis_current = Vector{Int64}
        if length(combinations) <= num_kopt_swaps_upperbound
            num_swaps = length(combinations)
        else
            num_swaps = num_kopt_swaps_upperbound
        end
        for i in 1:num_swaps
            if adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1])] == 0 && 
            adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1])] == 0
                G = add_multiple_edges!(G, [(src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1])),(src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1]))])
                if is_cyclic(G)
                    cycle_basis_current = cycle_basis(G)
                    algebraic_connectivity_tracker = 0.0
                    vertices_tracker = [(undef, undef), (undef, undef)]
                    algebraic_connectivity_tracker, vertices_tracker = vertices_tracker_span_update_two_edges!(G, adjacency_augment_graph, 0, 0, algebraic_connectivity_tracker, vertices_tracker)
                    G = rem_multiple_edges!(G, vertices_tracker)
                end
            end
        end
        updated_adjacency_graph_ac_tuple = update_kopt_adjacency!(G, adjacency_augment_graph, adjacency_graph_ac_tuple)
        
    elseif kopt_parameter == 3
        combinations = edge_combinations(length(sorted_edge_edgeweight_list), kopt_parameter)
        cycle_basis_current = Vector{Int64}
        if length(combinations) <= num_kopt_swaps_upperbound
            num_swaps = length(combinations)
        else
            num_swaps = num_kopt_swaps_upperbound
        end
        for i in 1:num_swaps
            if adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1])] == 0 &&  
            adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1])] == 0 &&  
            adjacency_matrix(G)[src(sorted_edge_edgeweight_list[combinations[i][3]][1]), dst(sorted_edge_edgeweight_list[combinations[i][3]][1])] == 0
                G = add_multiple_edges!(G, [(src(sorted_edge_edgeweight_list[combinations[i][1]][1]), dst(sorted_edge_edgeweight_list[combinations[i][1]][1])),(src(sorted_edge_edgeweight_list[combinations[i][2]][1]), dst(sorted_edge_edgeweight_list[combinations[i][2]][1])),(src(sorted_edge_edgeweight_list[combinations[i][3]][1]), dst(sorted_edge_edgeweight_list[combinations[i][3]][1]))])
                if is_cyclic(G)
                    algebraic_connectivity_tracker = 0.0
                    vertices_tracker = [(undef, undef), (undef, undef), (undef, undef)]
                    for j in 1:(length(cycle_basis(G)[3])-1)
                        for k in (j+1):length(cycle_basis(G)[3])
                            if adjacency_matrix(G)[cycle_basis(G)[3][j], cycle_basis(G)[3][k]] == 1
                                algebraic_connectivity_tracker, vertices_tracker = vertices_tracker_span_update_two_edges!(G, adjacency_augment_graph, j, k, algebraic_connectivity_tracker, vertices_tracker)
                            end
                        end
                    end
                    G = rem_multiple_edges!(G, vertices_tracker)
                end
            end
        end
        updated_adjacency_graph_ac_tuple = update_kopt_adjacency!(G, adjacency_augment_graph, adjacency_graph_ac_tuple)
    end
    return updated_adjacency_graph_ac_tuple
end

function vertices_tracker_span_update_two_edges!(
    G,
    adjacency_augment_graph,
    j, # Index for third edge src if kopt_parameter is 3
    k, # Index for third edge dst if kopt_parameter is 3
    algebraic_connectivity_tracker,
    vertices_tracker,
    )
    cycle_basis_current = cycle_basis(G)
    for l in 1:(length(cycle_basis_current[2])-1)
        for m in (l+1):length(cycle_basis_current[2])
            for p in 1:(length(cycle_basis_current[1])-1)
                for q in (p+1):length(cycle_basis_current[1])
                    cycle_basis_to_check = []
                    if j == 0 && k == 0
                        cycle_basis_to_check = [(cycle_basis_current[2][l], cycle_basis_current[2][m]),(cycle_basis_current[1][p], cycle_basis_current[1][q])]
                        cycle_basis_condition = (cycle_basis_current[1][p], cycle_basis_current[1][q]) !== (cycle_basis_current[2][l], cycle_basis_current[2][m])
                    else
                        cycle_basis_to_check = [(cycle_basis_current[3][j], cycle_basis_current[3][k]),(cycle_basis_current[2][l], cycle_basis_current[2][m]),(cycle_basis_current[1][p], cycle_basis_current[1][q])]
                        cycle_basis_condition = (cycle_basis_current[3][j], cycle_basis_current[3][k]) !== (cycle_basis_current[2][l], cycle_basis_current[2][m]) &&
                            (cycle_basis_current[3][j], cycle_basis_current[3][k]) !== (cycle_basis_current[1][p], cycle_basis_current[1][q]) && 
                            (cycle_basis_current[1][p], cycle_basis_current[1][q]) !== (cycle_basis_current[2][l], cycle_basis_current[2][m])
                    end
                    if cycle_basis_condition
                        if  adjacency_matrix(G)[cycle_basis_current[2][l], cycle_basis_current[2][m]] == 1 && adjacency_matrix(G)[cycle_basis_current[1][p], cycle_basis_current[1][q]] == 1
                            G = rem_multiple_edges!(G, cycle_basis_to_check)
                            if algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G))) > algebraic_connectivity_tracker
                                algebraic_connectivity_tracker = algebraic_connectivity(adjacency_augment_graph .* Matrix(adjacency_matrix(G)))
                                vertices_tracker = cycle_basis_to_check
                            end
                            G = add_multiple_edges!(G, cycle_basis_to_check)
                        end
                    end
                end
            end
        end
    end
    return algebraic_connectivity_tracker, vertices_tracker
end

function refinement_tree(
    G, 
    adjacency_base_graph,
    adjacency_augment_graph,
    sorted_edge_fiedlerweight_list, 
    kopt_parameter,
    num_kopt_swaps_upperbound,
    )

    if kopt_parameter == 1
        for i in eachindex(sorted_edge_fiedlerweight_list)
            if !(has_edge(G, src(sorted_edge_fiedlerweight_list[i][1]), dst(sorted_edge_fiedlerweight_list[i][1])))
                add_edge!(G, src(sorted_edge_fiedlerweight_list[i][1]), dst(sorted_edge_fiedlerweight_list[i][1]))
                algebraic_connectivity_tracker = 0.0
                vertices_tracker = (undef, undef)
                for edge in edges(G)
                    if adjacency_base_graph[src(edge), dst(edge)] == 0
                        rem_edge!(G, src(edge), dst(edge))
                        if  algebraic_connectivity((adjacency_augment_graph + adjacency_base_graph) .* Matrix(adjacency_matrix(G))) > algebraic_connectivity_tracker
                            algebraic_connectivity_tracker =  algebraic_connectivity((adjacency_augment_graph + adjacency_base_graph) .* Matrix(adjacency_matrix(G)))
                            vertices_tracker = edge
                        end
                        add_edge!(G, src(edge), dst(edge))
                    end
                end
                rem_edge!(G, vertices_tracker)
            end
        end

    elseif kopt_parameter == 2
        combinations = edge_combinations(length(sorted_edge_fiedlerweight_list), kopt_parameter)
        if length(combinations) <= num_kopt_swaps_upperbound
            num_swaps = length(combinations)
        else
            num_swaps = num_kopt_swaps_upperbound
        end
        for i in 1:num_swaps
            if !(has_edge(G, src(sorted_edge_fiedlerweight_list[combinations[i][1]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][1]][1])))  && 
            !(has_edge(G, src(sorted_edge_fiedlerweight_list[combinations[i][2]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][2]][1])))
                G = add_multiple_edges!(G, [(src(sorted_edge_fiedlerweight_list[combinations[i][1]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][1]][1])),(src(sorted_edge_fiedlerweight_list[combinations[i][2]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][2]][1]))])
                algebraic_connectivity_tracker = 0.0
                vertices_tracker = [(undef, undef),(undef, undef)]
                algebraic_connectivity_tracker, vertices_tracker = vertices_tracker_update_two_edges!(G, adjacency_base_graph, adjacency_augment_graph, 0, algebraic_connectivity_tracker,vertices_tracker)
                G = rem_multiple_edges!(G, vertices_tracker)
            end
        end

    elseif kopt_parameter == 3
        combinations = edge_combinations(length(sorted_edge_fiedlerweight_list), kopt_parameter)
        if length(combinations) <= num_kopt_swaps_upperbound
            num_swaps = length(combinations)
        else
            num_swaps = num_kopt_swaps_upperbound
        end
        for i in 1:num_swaps
            if !(has_edge(G, src(sorted_edge_fiedlerweight_list[combinations[i][1]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][1]][1])))  && 
            !(has_edge(G, src(sorted_edge_fiedlerweight_list[combinations[i][2]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][2]][1]))) &&
            !(has_edge(G, src(sorted_edge_fiedlerweight_list[combinations[i][3]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][3]][1]))) 
                G = add_multiple_edges!(G, [(src(sorted_edge_fiedlerweight_list[combinations[i][1]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][1]][1])),(src(sorted_edge_fiedlerweight_list[combinations[i][2]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][2]][1])), (src(sorted_edge_fiedlerweight_list[combinations[i][3]][1]), dst(sorted_edge_fiedlerweight_list[combinations[i][3]][1]))]) 
                algebraic_connectivity_tracker = 0.0
                vertices_tracker = [(undef, undef), (undef, undef), (undef, undef)]
                for j in eachindex(length(collect(edges(G))) - 2)
                    if (adjacency_base_graph[src(edge_list[j]), dst(edge_list[j])] == 0)
                        algebraic_connectivity_tracker, vertices_tracker = vertices_tracker_update_two_edges!(G, adjacency_base_graph, adjacency_augment_graph, j, algebraic_connectivity_tracker,vertices_tracker)
                    end
                end
                G = rem_multiple_edges!(G, vertices_tracker)
            end
        end
    end

    return Matrix(adjacency_matrix(G)), algebraic_connectivity((adjacency_augment_graph + adjacency_base_graph) .* Matrix(adjacency_matrix(G)))
end

function vertices_tracker_update_two_edges!(
    G, 
    adjacency_base_graph,
    adjacency_augment_graph,
    j, # Index for third edge if kopt_parameter is 3
    algebraic_connectivity_tracker, 
    vertices_tracker
    )
    edge_list = collect(edges(G))
    for k in j+1:(length(edge_list)-1)
        for l in k+1:length(edge_list)
            edges_to_check = []
            if j == 0
                edges_to_check = [(src(edge_list[k]), dst(edge_list[k])), (src(edge_list[l]), dst(edge_list[l]))]
            else
                edges_to_check = [(src(edge_list[j]), dst(edge_list[j])),(src(edge_list[k]), dst(edge_list[k])), (src(edge_list[l]), dst(edge_list[l]))]
            end
            if (adjacency_base_graph[src(edge_list[k]), dst(edge_list[k])] == 0) && (adjacency_base_graph[src(edge_list[l]), dst(edge_list[l])] == 0)
                G = rem_multiple_edges!(G, edges_to_check)
                if  algebraic_connectivity((adjacency_augment_graph + adjacency_base_graph) .* Matrix(adjacency_matrix(G))) > algebraic_connectivity_tracker
                    algebraic_connectivity_tracker =  algebraic_connectivity((adjacency_augment_graph + adjacency_base_graph) .* Matrix(adjacency_matrix(G)))
                    vertices_tracker = edges_to_check
                end
                G = add_multiple_edges!(G, edges_to_check)
            end
        end
    end
    
    return algebraic_connectivity_tracker, vertices_tracker
end

