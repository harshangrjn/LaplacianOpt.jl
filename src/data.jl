function get_data(params::Dict{String, Any})
    
    if "data_dict" in keys(params)
        data_dict = params["data_dict"]
    else 
        Memento.error(_LOGGER, "Data dictionary has to be specified in input params")
    end

    num_nodes = data_dict["num_nodes"] 
    instance  = params["instance"]

    if "augment_budget" in keys(params)
        if params["augment_budget"] <= data_dict["num_edges_to_augment"]
            augment_budget = params["augment_budget"]
        else
            augment_budget = data_dict["num_edges_to_augment"]
            Memento.info(_LOGGER, "Setting edge augmentation budget to $(data_dict["num_edges_to_augment"])")
        end
    else
        if data_dict["num_edges_to_augment"] == (num_nodes*(num_nodes-1))/2
            augment_budget = num_nodes-1
            Memento.info(_LOGGER, "Setting edge augmentation budget to a spanning tree")
        else
            augment_budget = data_dict["num_edges_to_augment"]
            Memento.info(_LOGGER, "Setting edge augmentation budget to $(data_dict["num_edges_to_augment"])")
        end
    end

    # Solution type 
    if "solution_type" in keys(params)
        solution_type = params["solution_type"]
    else
        solution_type = "exact"
    end

    # Relax Integrality 
    if "relax_integrality" in keys(params)
        relax_integrality = params["relax_integrality"]
    else
        # default value
        relax_integrality = false
    end

    # Tolerance to verify zero values
    if "tol_zero" in keys(params)
        tol_zero = params["tol_zero"]
    else
        # default value
        tol_zero = 1E-6
    end

    # Tolerance to verify PSD-ness of a matrix
    if "tol_psd" in keys(params)
        tol_psd = params["tol_psd"]
    else
        # default value
        tol_psd = 1E-6
    end

    if "eigen_cuts_full" in keys(params)
        eigen_cuts_full = params["eigen_cuts_full"]
        if eigen_cuts_full
            Memento.info(_LOGGER, "Applying full-sized eigen cuts")
        end
    else
        #default value
        Memento.info(_LOGGER, "Applying full-sized eigen cuts")
        eigen_cuts_full = true
    end

    if "soc_linearized_cuts" in keys(params)
        soc_linearized_cuts = params["soc_linearized_cuts"]
        if soc_linearized_cuts
            Memento.info(_LOGGER, "Applying linearized SOC cuts (2x2 minors)")
        end
    else
        #default value
        soc_linearized_cuts = false
    end

    if "topology_flow_cuts" in keys(params)
        topology_flow_cuts = params["topology_flow_cuts"]
    else
        #default value
        Memento.info(_LOGGER, "Applying topology flow cuts")
        topology_flow_cuts = true
    end

    if eigen_cuts_full || topology_flow_cuts || soc_linearized_cuts
        lazy_callback_status = true
    end

    if "lazycuts_logging" in keys(params)
        lazycuts_logging = params["lazycuts_logging"]
    else
        #default value
        lazycuts_logging = false
    end

    data = Dict{String, Any}("num_nodes"               => num_nodes,
                             "instance"                => instance,
                             "num_edges_existing"      => data_dict["num_edges_existing"],
                             "num_edges_to_augment"    => data_dict["num_edges_to_augment"],
                             "augment_budget"          => augment_budget,
                             "adjacency_base_graph"    => data_dict["adjacency_base_graph"],
                             "adjacency_augment_graph" => data_dict["adjacency_augment_graph"],
                             "solution_type"           => solution_type,
                             "tol_zero"                => tol_zero,
                             "tol_psd"                 => tol_psd,
                             "eigen_cuts_full"         => eigen_cuts_full,
                             "soc_linearized_cuts"     => soc_linearized_cuts,
                             "lazycuts_logging"        => lazycuts_logging,
                             "topology_flow_cuts"      => topology_flow_cuts,
                             "lazy_callback_status"    => lazy_callback_status,
                             "relax_integrality"       => relax_integrality)

    # Optimizer
    if "optimizer" in keys(params)
        data["optimizer"] = params["optimizer"]
    end

    return data
end

function parse_file(file_path::String)
    data_dict = JSON.parsefile(file_path)
    
    if "num_nodes" in keys(data_dict)
        num_nodes = data_dict["num_nodes"]
    else
        Memento.error(_LOGGER, "Number of nodes is missing in the input data file")
    end

    adjacency_base_graph    = Matrix{Float64}(zeros(num_nodes, num_nodes))
    adjacency_augment_graph = Matrix{Float64}(zeros(num_nodes, num_nodes))
    
    num_edges_existing   = 0 
    num_edges_to_augment = 0 

    for k=1:length(data_dict["edges_existing"])
        i, j = data_dict["edges_existing"][k][1]
        w_ij = data_dict["edges_existing"][k][2]

        if isapprox(abs(w_ij), 0, atol = 1E-6)
            w_ij = 0
        end

        LOpt._catch_data_input_error(num_nodes, i, j, w_ij)
        
        adjacency_base_graph[i,j] = w_ij
        adjacency_base_graph[j,i] = w_ij
        
        num_edges_existing += 1
    end

    for k=1:length(data_dict["edges_to_augment"])
        i, j = data_dict["edges_to_augment"][k][1]
        w_ij = data_dict["edges_to_augment"][k][2]

        if isapprox(abs(w_ij), 0, atol = 1E-6)
            w_ij = 0
        end

        LOpt._catch_data_input_error(num_nodes, i, j, w_ij)
        
        adjacency_augment_graph[i,j] = w_ij
        adjacency_augment_graph[j,i] = w_ij
        num_edges_to_augment += 1
    end

    if num_edges_to_augment == 0
        Memento.error(_LOGGER, "LaplacianOpt needs at least one edge to be able to augment to the base graph")
    # Assuming undirected graph
    elseif num_edges_existing == num_nodes*(num_nodes-1)/2 
        Memento.error(_LOGGER, "Input graph is already a complete graph; augmentation is unnecessary")
    end

    # Check for mutually exclusive sets of edges between base and augment graph
    if maximum((adjacency_augment_graph .> 0) + (adjacency_base_graph .> 0)) > 1
        Memento.error(_LOGGER, "Edge sets of base and augment graphs have to be mutually exclusive")
    end

    data_dict_new = Dict{String, Any}("num_nodes"               => num_nodes,
                                      "num_edges_existing"      => num_edges_existing,
                                      "num_edges_to_augment"    => num_edges_to_augment,
                                      "adjacency_base_graph"    => adjacency_base_graph,
                                      "adjacency_augment_graph" => adjacency_augment_graph)
    
    return data_dict_new
end

function _catch_data_input_error(num_nodes::Int64, i::Int64, j::Int64, w_ij::Number)
    if (i > num_nodes) || (j > num_nodes)
        Memento.error(_LOGGER, "Node pair ($i,$j) does not match with total number of nodes, $num_nodes")
    end

    if !(isapprox(abs(w_ij), 0, atol=1E-6)) & (w_ij < 0)
        Memento.error(_LOGGER, "LaplacianOpt does not support graphs with negative weights")
    end
end