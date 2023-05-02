"""
    get_data(params::Dict{String, Any})

Given the input params, this function preprocesses the data, catches any error and missing data, 
and outputs a detailed dictionary which forms the basis for building an optimization model. 
"""
function get_data(params::Dict{String,Any})
    _header_color = :cyan
    printstyled(
        "LaplacianOpt version: ",
        Pkg.TOML.parse(read(string(pkgdir(LOpt), "/Project.toml"), String))["version"],
        "\n";
        color = _header_color,
        bold = true,
    )

    if "data_dict" in keys(params)
        data_dict = params["data_dict"]
    else
        Memento.error(_LOGGER, "Data dictionary has to be specified in input params")
    end

    data_dict = LOpt._pre_process_data(data_dict)

    Memento.info(
        _LOGGER,
        "Number of nodes in the graph: $(data_dict["num_nodes"])",
    )

    if "augment_budget" in keys(params)
        if params["augment_budget"] <= data_dict["num_edges_to_augment"]
            augment_budget = params["augment_budget"]
        else
            augment_budget = data_dict["num_edges_to_augment"]
            Memento.info(
                _LOGGER,
                "Setting edge augmentation budget to $(data_dict["num_edges_to_augment"])",
            )
        end
    else
        augment_budget = data_dict["num_edges_to_augment"]
        Memento.info(
            _LOGGER,
            "Setting edge augmentation budget to $(data_dict["num_edges_to_augment"])",
        )
    end

    Memento.info(
        _LOGGER,
        "Number of edges in the base graph: $(data_dict["num_edges_existing"])",
    )
    Memento.info(
        _LOGGER,
        "Number of candidate edges to augment: $(data_dict["num_edges_to_augment"])",
    )
    Memento.info(_LOGGER, "Augment budget: $(augment_budget)")

    # Graph type
    if "graph_type" in keys(params)
        graph_type = params["graph_type"]
        if !(graph_type in ["any", "hamiltonian_cycle"])
            Memento.error(_LOGGER, "Invalid graph type, $graph_type, in input params")
        end
    else
        # default value
        graph_type = "any"
    end

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

    if "num_nodes" in keys(data_dict)
        num_nodes = data_dict["num_nodes"]
    else
        Memento.error(_LOGGER, "Number of nodes is missing in the input data file")
    end

    adjacency_base_graph = Matrix{Float64}(zeros(num_nodes, num_nodes))
    adjacency_augment_graph = Matrix{Float64}(zeros(num_nodes, num_nodes))

    for k in 1:length(data_dict["edges_existing"])
        i, j = data_dict["edges_existing"][k][1]
        w_ij = data_dict["edges_existing"][k][2]

        if isapprox(abs(w_ij), 0, atol = 1E-6)
            w_ij = 0
        end

        LOpt._catch_data_input_error(num_nodes, i, j, w_ij)

        adjacency_base_graph[i, j] = w_ij
        adjacency_base_graph[j, i] = w_ij
    end

    for k in 1:length(data_dict["edges_to_augment"])
        i, j = data_dict["edges_to_augment"][k][1]
        w_ij = data_dict["edges_to_augment"][k][2]

        if isapprox(abs(w_ij), 0, atol = 1E-6)
            w_ij = 0
        end

        LOpt._catch_data_input_error(num_nodes, i, j, w_ij)

        adjacency_augment_graph[i, j] = w_ij
        adjacency_augment_graph[j, i] = w_ij
    end

    data_dict_new = Dict{String,Any}(
        "num_nodes" => num_nodes,
        "adjacency_base_graph" => adjacency_base_graph,
        "adjacency_augment_graph" => adjacency_augment_graph,
    )

    return data_dict_new
end

function _pre_process_data(data_dict::Dict{String,Any})
    num_edges_existing = 0
    num_edges_to_augment = 0

    num_nodes = data_dict["num_nodes"]
    adjacency_base_graph = data_dict["adjacency_base_graph"]
    adjacency_augment_graph = data_dict["adjacency_augment_graph"]

    for i in 1:(num_nodes-1)
        for j in (i+1):num_nodes
            if isapprox(abs(adjacency_base_graph[i, j]), 0, atol = 1E-6)
                adjacency_base_graph[i, j] = 0
            end
            if isapprox(abs(adjacency_augment_graph[i, j]), 0, atol = 1E-6)
                adjacency_augment_graph[i, j] = 0
            end
            if !(isapprox(
                adjacency_base_graph[i, j],
                adjacency_base_graph[j, i],
                atol = 1E-5,
            ))
                Memento.error("Adjacency matrix of the base graph has to be symmetric")
            end
            if !(isapprox(
                adjacency_augment_graph[i, j],
                adjacency_augment_graph[j, i],
                atol = 1E-5,
            ))
                Memento.error("Adjacency matrix of the augment graph has to be symmetric")
            end
            if !(isapprox(adjacency_base_graph[i, j], 0, atol = 1E-6))
                num_edges_existing += 1
            end
            if !(isapprox(adjacency_augment_graph[i, j], 0, atol = 1E-6))
                num_edges_to_augment += 1
            end
        end
    end

    # Base graph connectivity
    is_base_graph_connected = false
    if num_edges_existing > 0
        !(isapprox(
            abs(LOpt.algebraic_connectivity(adjacency_base_graph)),
            0,
            atol = 1E-6,
        )) && (is_base_graph_connected = true)
    end

    data_dict_new = Dict{String,Any}()
    data_dict_new["num_nodes"] = num_nodes
    data_dict_new["adjacency_base_graph"] = adjacency_base_graph
    data_dict_new["adjacency_augment_graph"] = adjacency_augment_graph
    data_dict_new["num_edges_existing"] = num_edges_existing
    data_dict_new["num_edges_to_augment"] = num_edges_to_augment
    data_dict_new["is_base_graph_connected"] = is_base_graph_connected

    return data_dict_new
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

    if data["num_edges_to_augment"] == 0
        Memento.error(
            _LOGGER,
            "At least one edge has to be specified to be able to augment to the base graph",
        )
        # Assuming undirected graph
    elseif data["num_edges_existing"] == num_nodes * (num_nodes - 1) / 2
        Memento.error(
            _LOGGER,
            "Input graph is already a complete graph; augmentation is unnecessary",
        )
    elseif (data["num_edges_existing"] == 0) && (data["augment_budget"] < (num_nodes - 1))
        Memento.error(
            _LOGGER,
            "Detected trivial solutions with disconnected graphs. `augment_budget` may be insufficient.",
        )
    elseif !isinteger(data["augment_budget"]) || (data["augment_budget"] < -1E-6)
        Memento.error(_LOGGER, "Edge augmentation budget has to be a positive integer")
    end

    # Detect mutually exclusive sets of edges between base and augment graph
    if maximum(
        (data["adjacency_augment_graph"] .> 0) + (data["adjacency_base_graph"] .> 0),
    ) > 1
        Memento.error(
            _LOGGER,
            "Edge sets of base and augment graphs have to be mutually exclusive",
        )
    end

    # Detect free vertices 
    A = data["adjacency_augment_graph"] + data["adjacency_base_graph"]
    for i in 1:num_nodes
        if isapprox(sum(A[i, :]), 0, atol = 1E-6)
            Memento.error(
                _LOGGER,
                "Detected trivial solutions with disconnected graphs due to free vertices.",
            )
        end
    end

    # Detect tour infeasibility
    if (data["graph_type"] == "hamiltonian_cycle") && (data["num_edges_existing"] == 0)
        if !(data["num_edges_to_augment"] >= data["num_nodes"]) ||
           !(data["augment_budget"] == data["num_nodes"])
            Memento.error(
                _LOGGER,
                "Detected infeasibility due to the number of augmentation edges incompatible for a hamiltonian cycle",
            )
        end
    end
end
