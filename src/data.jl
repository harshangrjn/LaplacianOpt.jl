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

    # Tolerance to verify zero values
    if "tol_zero" in keys(params)
        tol_zero = params["tol_zero"]

        if params["tol_zero"] >= 1E-4
            Memento.warn(_LOGGER, "Zero tolerance value may be too high for numerical accuracy in solutions")
        end

    else
        # default value
        tol_zero = 1E-6
    end

    # Tolerance to verify PSD-ness of a matrix
    if "tol_psd" in keys(params)
        tol_psd = params["tol_psd"]

        if params["tol_psd"] >= 1E-3
            Memento.warn(_LOGGER, "PSD tolerance value may be too high for PSD feasibility")
        end
        
    else
        # default value
        tol_psd = 1E-6
    end

    if "eigen_cuts_full" in keys(params)
        eigen_cuts_full = params["eigen_cuts_full"]
    else
        #default value
        Memento.info(_LOGGER, "Turning on full-sized eigen cuts")
        eigen_cuts_full = true
    end

    if "topology_flow_cuts" in keys(params)
        topology_flow_cuts = params["topology_flow_cuts"]
    else
        #default value
        topology_flow_cuts = true
    end

    if eigen_cuts_full || topology_flow_cuts
        lazy_callback_status = true
    end

    if "lazycuts_logging" in keys(params)
        lazycuts_logging = params["lazycuts_logging"]
    else
        #default value
        lazycuts_logging = false
    end

    data = Dict{String, Any}("num_nodes" => n,
                             "instance" => instance,
                             "edge_weights" => edge_weigths_matrix,
                             "tol_zero" => tol_zero,
                             "tol_psd" => tol_psd,
                             "eigen_cuts_full" => eigen_cuts_full,
                             "lazycuts_logging" => lazycuts_logging,
                             "topology_flow_cuts" => topology_flow_cuts,
                             "lazy_callback_status" => lazy_callback_status,
                             "relax_integrality" => params["relax_integrality"])

    return data
end

function convert_array_to_matrix(n::Int, edge_weights::Any)
    
    instance_matrix = Matrix{Float64}(zeros(n,n))

    for k = 1:length(edge_weights)
        i, j = edge_weights[k][1]
        w_ij = edge_weights[k][2]

        if w_ij < 0 
            Memento.error(_LOGGER, "LaplacianOpt does not support graphs with negative weights")
        end

        instance_matrix[i,j] = w_ij
        instance_matrix[j,i] = w_ij
    end

    return instance_matrix
end