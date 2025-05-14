"""
    get_data(params::Dict{String, Any})

Given the input params, this function preprocesses the data, catches any error and missing data, 
and outputs a detailed dictionary which forms the basis for building an optimization model. 
"""
function get_data(params::Dict{String,Any})
    # version info
    _header_color = :cyan
    printstyled(
        "LaplacianOpt version: ",
        Pkg.TOML.parse(read(string(pkgdir(LOpt), "/Project.toml"), String))["version"],
        "\n";
        color = _header_color,
        bold = true,
    )

    # data dictionary
    if !haskey(params, "data_dict")
        Memento.error(_LOGGER, "Data dictionary has to be specified in input params")
    end
    data_dict = LOpt._pre_process_data(params["data_dict"])

    Memento.info(_LOGGER, "Number of nodes in the graph: $(data_dict["num_nodes"])")

    # augment budget
    augment_budget =
        if haskey(params, "augment_budget") &&
           params["augment_budget"] <= data_dict["num_edges_to_augment"]
            params["augment_budget"]
        else
            budget = data_dict["num_edges_to_augment"]
            Memento.info(_LOGGER, "Setting edge augmentation budget to $budget")
            budget
        end

    # Log graph information
    Memento.info(
        _LOGGER,
        "Number of edges in the base graph: $(data_dict["num_edges_existing"])",
    )
    Memento.info(
        _LOGGER,
        "Number of candidate edges to augment: $(data_dict["num_edges_to_augment"])",
    )
    Memento.info(_LOGGER, "Augment budget: $(augment_budget)")

    # graph type
    graph_type = "any"  # default
    if haskey(params, "graph_type")
        if params["graph_type"] in ["any", "hamiltonian_cycle"]
            graph_type = params["graph_type"]
        else
            Memento.error(
                _LOGGER,
                "Invalid graph type, $(params["graph_type"]), in input params",
            )
        end
    end

    # data dictionary
    data = Dict{String,Any}(
        "num_nodes" => data_dict["num_nodes"],
        "num_edges_existing" => data_dict["num_edges_existing"],
        "num_edges_to_augment" => data_dict["num_edges_to_augment"],
        "augment_budget" => augment_budget,
        "adjacency_base_graph" => data_dict["adjacency_base_graph"],
        "adjacency_augment_graph" => data_dict["adjacency_augment_graph"],
        "is_base_graph_connected" => data_dict["is_base_graph_connected"],
        "graph_type" => graph_type,
    )

    LOpt._detect_infeasbility_in_data(data)
    return data
end

function parse_file(file_path::String)
    data_dict = JSON.parsefile(file_path)

    haskey(data_dict, "num_nodes") ||
        Memento.error(_LOGGER, "Number of nodes is missing in the input data file")
    num_nodes = data_dict["num_nodes"]

    adjacency_base_graph = zeros(Float64, num_nodes, num_nodes)
    adjacency_augment_graph = zeros(Float64, num_nodes, num_nodes)

    # Process existing edges
    for edge in data_dict["edges_existing"]
        (i, j), w_ij = edge
        w_ij = isapprox(abs(w_ij), 0, atol = 1E-6) ? 0 : w_ij

        LOpt._catch_data_input_error(num_nodes, i, j, w_ij)
        adjacency_base_graph[i, j] = adjacency_base_graph[j, i] = w_ij
    end

    # Process augment edges
    for edge in data_dict["edges_to_augment"]
        (i, j), w_ij = edge
        w_ij = isapprox(abs(w_ij), 0, atol = 1E-6) ? 0 : w_ij

        LOpt._catch_data_input_error(num_nodes, i, j, w_ij)
        adjacency_augment_graph[i, j] = adjacency_augment_graph[j, i] = w_ij
    end

    return Dict{String,Any}(
        "num_nodes" => num_nodes,
        "adjacency_base_graph" => adjacency_base_graph,
        "adjacency_augment_graph" => adjacency_augment_graph,
    )
end

function _pre_process_data(data_dict::Dict{String,Any})
    num_nodes = data_dict["num_nodes"]
    adjacency_base_graph = data_dict["adjacency_base_graph"]
    adjacency_augment_graph = data_dict["adjacency_augment_graph"]

    num_edges_existing = 0
    num_edges_to_augment = 0

    # Process matrices in one pass
    for i in 1:(num_nodes-1), j in (i+1):num_nodes
        if isapprox(abs(adjacency_base_graph[i, j]), 0, atol = 1E-6)
            adjacency_base_graph[i, j] = adjacency_base_graph[j, i] = 0
        elseif !isapprox(
            adjacency_base_graph[i, j],
            adjacency_base_graph[j, i],
            atol = 1E-5,
        )
            Memento.error("Adjacency matrix of the base graph has to be symmetric")
        else
            num_edges_existing += 1
        end

        if isapprox(abs(adjacency_augment_graph[i, j]), 0, atol = 1E-6)
            adjacency_augment_graph[i, j] = adjacency_augment_graph[j, i] = 0
        elseif !isapprox(
            adjacency_augment_graph[i, j],
            adjacency_augment_graph[j, i],
            atol = 1E-5,
        )
            Memento.error("Adjacency matrix of the augment graph has to be symmetric")
        else
            num_edges_to_augment += 1
        end
    end

    # Check base graph connectivity
    is_base_graph_connected =
        num_edges_existing > 0 &&
        !isapprox(abs(LOpt.algebraic_connectivity(adjacency_base_graph)), 0, atol = 1E-6)

    return Dict{String,Any}(
        "num_nodes" => num_nodes,
        "adjacency_base_graph" => adjacency_base_graph,
        "adjacency_augment_graph" => adjacency_augment_graph,
        "num_edges_existing" => num_edges_existing,
        "num_edges_to_augment" => num_edges_to_augment,
        "is_base_graph_connected" => is_base_graph_connected,
    )
end

"""
    _catch_data_input_error(num_nodes::Int64, i::Int64, j::Int64, w_ij::Number)

Given the number of nodes, from and to nodes of an edge, and an edge weight, this function catches 
any input data error in the JSON file. 
"""
function _catch_data_input_error(num_nodes::Int64, i::Int64, j::Int64, w_ij::Number)
    if (i > num_nodes) || (j > num_nodes)
        Memento.error(
            _LOGGER,
            "Node pair ($i,$j) does not match with total number of nodes, $num_nodes",
        )
    end

    if !(isapprox(abs(w_ij), 0, atol = 1E-6)) & (w_ij < 0)
        Memento.error(_LOGGER, "Graphs with negative weights are not supported")
    end
end

"""
    _detect_infeasbility_in_data(data::Dict{String, Any})

Given the pre-processed data dictionary, this function detects any infeasibility before 
building the optimization model. 
"""
function _detect_infeasbility_in_data(data::Dict{String,Any})
    num_nodes = data["num_nodes"]
    num_edges_existing = data["num_edges_existing"]
    num_edges_to_augment = data["num_edges_to_augment"]
    augment_budget = data["augment_budget"]
    graph_type = get(data, "graph_type", "")

    # Check for basic infeasibility conditions
    num_edges_to_augment == 0 &&
        Memento.error(_LOGGER, "At least one edge must be available for augmentation")
    num_edges_existing == num_nodes * (num_nodes - 1) / 2 && Memento.error(
        _LOGGER,
        "Input graph is already complete; augmentation unnecessary",
    )
    (num_edges_existing == 0 && augment_budget < (num_nodes - 1)) &&
        Memento.error(_LOGGER, "Insufficient augment_budget for connectivity")
    (!isinteger(augment_budget) || augment_budget < -1E-6) &&
        Memento.error(_LOGGER, "Edge augmentation budget must be positive integer")

    # Check for mutually exclusive edge sets
    maximum(
        (data["adjacency_augment_graph"] .> 0) + (data["adjacency_base_graph"] .> 0),
    ) > 1 && Memento.error(
        _LOGGER,
        "Edge sets of base and augment graphs must be mutually exclusive",
    )

    # Check for free vertices
    A = data["adjacency_augment_graph"] + data["adjacency_base_graph"]
    any(isapprox.(sum(A, dims = 2), 0, atol = 1E-6)) &&
        Memento.error(_LOGGER, "Detected free vertices resulting in disconnected graphs")

    # Check Hamiltonian cycle feasibility
    if graph_type == "hamiltonian_cycle" && num_edges_existing == 0
        (num_edges_to_augment < num_nodes || augment_budget != num_nodes) &&
            Memento.error(_LOGGER, "Insufficient edges for Hamiltonian cycle")
    end
end
