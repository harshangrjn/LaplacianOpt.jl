function get_data(params::Dict{String, Any})
    n = params["num_nodes"] 
    instance = params["instance"]

    file_path = joinpath(dirname(pathof(LaplacianOpt)),"..", "examples/instances/$(n)_nodes/$(n)_$(instance).json")
    
    data_dict = JSON.parsefile(file_path)

    if data_dict["num_nodes"] != n 
        Memento.error(_LOGGER, "Mismatch in number of input nodes")
    end

    if params["data_type"] == "old"
        edge_weights = data_dict["old_data"]["edge_weights"]
    elseif params["data_type"] == "new"
        edge_weights = data_dict["new_data"]["edge_weights"]
    end

    edge_weigths_matrix = convert_array_to_matrix(n, edge_weights)

    data = Dict{String, Any}("num_nodes" => n,
                             "instance" => instance,
                             "edge_weights" => edge_weigths_matrix,
                             "relax_integrality" => params["relax_integrality"])

    return data
end

function convert_array_to_matrix(n::Int, edge_weights::Any)
    
    instance_matrix = Matrix{Float64}(zeros(n,n))

    for k = 1:length(edge_weights)
        i, j = edge_weights[k][1]
        w_ij = edge_weights[k][2]

        if w_ij < 0 
            Memento.error(_LOGGER, "LaplaciaOpt does not support graphs with negative weights")
        end

        instance_matrix[i,j] = w_ij
        instance_matrix[j,i] = w_ij
    end

    return instance_matrix
end