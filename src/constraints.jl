#------------------------------------------------------#
# Build all constraints of the LaplacianOptModel here  #
#------------------------------------------------------#

function constraint_build_W_var_matrix(lom::LaplacianOptModel)
    num_nodes               = lom.data["num_nodes"]
    num_edges_existing      = lom.data["num_edges_existing"]
    adjacency_base_graph    = lom.data["adjacency_base_graph"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    adjacency_full_graph = adjacency_augment_graph
    (num_edges_existing > 0) && (adjacency_full_graph += adjacency_base_graph)
        
    # Diagonal entries
    JuMP.@constraint(lom.model, [i=1:num_nodes], lom.variables[:W_var][i,i] == sum(adjacency_full_graph[i,:] .* lom.variables[:z_var][i,:]) - lom.variables[:γ_var]*(num_nodes-1)/num_nodes)

    # Off-diagonal entries
    JuMP.@constraint(lom.model, [i=1:(num_nodes-1), j=(i+1):num_nodes], lom.variables[:W_var][i,j] == -adjacency_full_graph[i,j] * lom.variables[:z_var][i,j] + lom.variables[:γ_var]/num_nodes)

    return
end

function constraint_single_vertex_cutset(lom::LaplacianOptModel)
    num_nodes               = lom.data["num_nodes"]
    num_edges_existing      = lom.data["num_edges_existing"]
    adjacency_base_graph    = lom.data["adjacency_base_graph"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    adjacency_full_graph = adjacency_augment_graph
    (num_edges_existing > 0) && (adjacency_full_graph += adjacency_base_graph)

    for i = 1:num_nodes
        if sum(adjacency_full_graph[i,:]) > 0
            JuMP.@constraint(lom.model, sum(lom.variables[:z_var][i,j] for j in setdiff((1:num_nodes), i)) >= 1)
        end
    end
    
    return
end

function constraint_augment_edges_budget(lom::LaplacianOptModel)
    num_nodes               = lom.data["num_nodes"]
    augment_budget          = lom.data["augment_budget"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    JuMP.@constraint(lom.model, sum(lom.variables[:z_var][i,j] for i=1:(num_nodes-1), j=(i+1):num_nodes 
                                    if adjacency_augment_graph[i,j] > 0) == augment_budget)
    
    return
end

function constraint_lazycallback_wrapper(lom::LaplacianOptModel)
    
    function constraint_lazycallback_cuts(cb_cuts)

        status = JuMP.callback_node_status(cb_cuts, lom.model)

        # if status == MOI.CALLBACK_NODE_STATUS_FRACTIONAL
            
        #     if lom.data["optimizer"] != "glpk"
        #         Memento.info(_LOGGER, "Callback status: solution is fractional")
        #     end

        if status == MOI.CALLBACK_NODE_STATUS_UNKNOWN
            Memento.error(_LOGGER, "Callback status: unknown - solution may not be integer feasible")

        else status == MOI.CALLBACK_NODE_STATUS_INTEGER
            
            W_val = JuMP.callback_value.(Ref(cb_cuts), lom.variables[:W_var])

            if lom.data["eigen_cuts_full"]    
                LOpt.constraint_eigen_cuts_on_full_matrix(W_val, cb_cuts, lom)
            end

            if lom.data["soc_linearized_cuts"]
                LOpt.constraint_soc_cuts_on_2minors(W_val, cb_cuts, lom)
            end

            if lom.data["eigen_cuts_2minors"]
                LOpt.constraint_eigen_cuts_on_2minors(W_val, cb_cuts, lom)
            end

            if lom.data["eigen_cuts_3minors"]
                LOpt.constraint_eigen_cuts_on_3minors(W_val, cb_cuts, lom)
            end

            if !(lom.data["is_base_graph_connected"]) && (lom.data["topology_flow_cuts"])
                z_val = abs.(JuMP.callback_value.(Ref(cb_cuts), lom.variables[:z_var]))         
                LOpt.get_rounded_zeros_and_ones!(z_val, lom.data["tol_zero"])
                LOpt.constraint_topology_flow_cuts(z_val, cb_cuts, lom)
            end
        end 
    end

    if lom.data["lazy_callback_status"]
        MOI.set(lom.model, MOI.LazyConstraintCallback(), constraint_lazycallback_cuts)
    end

end

function constraint_soc_cuts_on_2minors(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)
    
    num_nodes = lom.data["num_nodes"]

    for i = 1:(num_nodes-1)
        for j = (i+1):num_nodes
            W_val_minor = W_val[i,j]^2 - W_val[i,i] * W_val[j,j]
            if  W_val_minor >= lom.data["tol_zero"]                
                coefficient_vec = [-(W_val[i,j])^2
                                   2 * W_val[i,j] * W_val[i,i]]
                
                variable_vec = [lom.variables[:W_var][i,i], lom.variables[:W_var][i,j]]
                
                # Taylor's approximation of rotated SOC constraint (x^2 <= y*z)
                if !isapprox(W_val[i,i], 0, atol=1E-6)
                    con = JuMP.@build_constraint(coefficient_vec[1]*variable_vec[1] + coefficient_vec[2]*variable_vec[2] 
                                                    <= W_val[i,i]^2 * lom.variables[:W_var][j,j])
                    
                    MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)   
                end
                
                if lom.data["lazycuts_logging"]
                    Memento.info(_LOGGER, "Polyhedral relaxation cuts: 2x2 minor linearized SOC")
                end    
            end
        end
    end

    return
end

function constraint_eigen_cuts_on_full_matrix(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)
    
    LOpt._add_eigen_cut_lazy(W_val, cb_cuts, lom, collect(1:lom.data["num_nodes"]))
    return
end

function constraint_eigen_cuts_on_2minors(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)
    
    num_nodes = lom.data["num_nodes"]
    
    if num_nodes >= 3
        for i = 1:(num_nodes-1)
            for j = (i+1):num_nodes
                LOpt._add_eigen_cut_lazy(W_val, cb_cuts, lom, [i,j])     
            end
        end
    end
    return
end

function constraint_eigen_cuts_on_3minors(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)
    
    num_nodes = lom.data["num_nodes"]

    if num_nodes >= 4
        for i = 1:(num_nodes-2)
            for j = (i+1):num_nodes
                for k = (j+1):num_nodes
                    LOpt._add_eigen_cut_lazy(W_val, cb_cuts, lom, [i,j,k])
                end     
            end
        end
    end
    return
end

function _add_eigen_cut_lazy(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel, idx::Vector{Int64})
    W_val_minor = W_val[idx, idx]
    W_var_minor = lom.variables[:W_var][idx, idx]

    violated_eigen_vec = LOpt._violated_eigen_vector(W_val_minor, tol = lom.data["tol_psd"])
    if violated_eigen_vec !== nothing 
        con = JuMP.@build_constraint(violated_eigen_vec' * W_var_minor * violated_eigen_vec >= 0)
        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)

        if lom.data["lazycuts_logging"]
            if length(idx) == 2
                Memento.info(_LOGGER, "Polyhedral relaxation cuts: 2x2 minor (eigen)")
            elseif length(idx) == 3
                Memento.info(_LOGGER, "Polyhedral relaxation cuts: 3x3 minor (eigen)")
            elseif length(idx) == lom.data["num_nodes"]
                Memento.info(_LOGGER, "Polyhedral relaxation cuts: nxn matrix (eigen)")
            end
        end
    end

end

function constraint_topology_flow_cuts(z_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)

    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    cc_lazy = Graphs.connected_components(Graphs.SimpleGraph(abs.(z_val)))
    
    max_cc  = 5 # increase this to any greater integer value and it works
    if length(cc_lazy) > max_cc
        Memento.info(_LOGGER, "Polyhedral relaxation: flow cuts not added for integer solutions with $(length(cc_lazy)) connected components")
    end

    if length(cc_lazy) in 2:max_cc
       
        for k = 1:(length(cc_lazy)-1)
            # From cutset  
            cutset_f = cc_lazy[k]        
            # To cutset
            ii = setdiff(1:length(cc_lazy), k)
            cutset_t = reduce(vcat, cc_lazy[ii])
            
            if LOpt._is_flow_cut_valid(cutset_f, cutset_t, adjacency_augment_graph)
                con = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cutset_f, j in cutset_t 
                                            if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)

                MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)
            end
        end

        if lom.data["lazycuts_logging"]
            Memento.info(_LOGGER, "Polyhedral relaxation cuts: graph connectivity (#cc = $(length(cc_lazy)))")
        end 
    end

    return 
end

function constraint_topology_multi_commodity_flow(lom::LaplacianOptModel)

    num_nodes = lom.data["num_nodes"]
    
    # Flow balance on source node 
    JuMP.@constraint(lom.model, [k=1:(num_nodes-1)], sum(lom.variables[:flow_var][1,2:num_nodes,k]) - sum(lom.variables[:flow_var][2:num_nodes,1,k]) == 1)
    # Flow balance on target node 
    JuMP.@constraint(lom.model, [k=1:(num_nodes-1)], sum(lom.variables[:flow_var][:,k+1,k]) - sum(lom.variables[:flow_var][k+1,:,k]) == 1) #(k-th commodity corresponds to k+1th node)
    # Flow balance on remaining nodes 
    JuMP.@constraint(lom.model, [k=1:(num_nodes-1), v=(2:num_nodes); k!=v-1], sum(lom.variables[:flow_var][:,v,k]) - sum(lom.variables[:flow_var][v,:,k]) == 0)
    # Flow <= capacity 
    JuMP.@constraint(lom.model, [k=1:(num_nodes-1), i=1:num_nodes, j=1:num_nodes, kdash=1:(num_nodes-1)], lom.variables[:flow_var][i,j,k] + lom.variables[:flow_var][j,i,kdash]' <= lom.variables[:z_var][i,j])
    
    return
end