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

function constraint_topology_vertex_cutset(lom::LaplacianOptModel)
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

            if lom.data["eigen_cuts_full"]
                W_val = JuMP.callback_value.(Ref(cb_cuts), lom.variables[:W_var])

                LOpt.constraint_eigen_cuts_on_full_matrix(W_val, cb_cuts, lom)
            end

            if lom.data["soc_linearized_cuts"]
                W_val = JuMP.callback_value.(Ref(cb_cuts), lom.variables[:W_var])

                LOpt.constraint_soc_linearized_cuts_on_2minors(W_val, cb_cuts, lom)
            end

            if lom.data["topology_flow_cuts"]
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

function constraint_eigen_cuts_on_full_matrix(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)
    
    if typeof(LA.eigvals(W_val)) != Vector{Float64}
        Memento.error(_LOGGER, "PSD matrix (W) cannot have complex eigenvalues")
    end

    if LA.eigmin(W_val) <= -lom.data["tol_psd"]
        
        if LA.eigvals(W_val)[1] == LA.eigmin(W_val)
            
            violated_eigen_vec = LA.eigvecs(W_val)[:,1]
            
            LOpt.get_rounded_zeros_and_ones!(violated_eigen_vec, lom.data["tol_zero"])

            con = JuMP.@build_constraint(violated_eigen_vec' * lom.variables[:W_var] * violated_eigen_vec >= 0)
            
            MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con) 
        
            if lom.data["lazycuts_logging"]
                Memento.info(_LOGGER, "Polyhedral relaxation cuts: full-sized eigen")
            end    

        else

            Memento.warn(_LOGGER, "Eigen cut corresponding to the negative eigenvalue could not be added")
        end
    end

    return
end

function constraint_soc_linearized_cuts_on_2minors(W_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)
    
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

function constraint_topology_flow_cuts(z_val::Matrix{Float64}, cb_cuts, lom::LaplacianOptModel)

    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    cc_lazy = Graphs.connected_components(Graphs.SimpleGraph(abs.(z_val)))

    if length(cc_lazy) == 2 
        
        con = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cc_lazy[1], j in cc_lazy[2] 
                                         if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)  

        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con)      
        
        if lom.data["lazycuts_logging"]
            Memento.info(_LOGGER, "Polyhedral relaxation cuts: graph connectivity (#cc = $(length(cc_lazy)))")
        end

    elseif length(cc_lazy) == 3

        con1 = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cc_lazy[1], j in union(cc_lazy[2], cc_lazy[3])
                                      if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)  
                                              
        con2 = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cc_lazy[2], j in union(cc_lazy[3], cc_lazy[1]) 
                                      if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)  

        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con1)       
        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con2)
        
        if lom.data["lazycuts_logging"]
            Memento.info(_LOGGER, "Polyhedral relaxation cuts: graph connectivity (#cc = $(length(cc_lazy)))")
        end    

    elseif length(cc_lazy) == 4

        con1 = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cc_lazy[1], j in union(cc_lazy[2], cc_lazy[3], cc_lazy[4])
                                      if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)  
                
        con2 = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cc_lazy[2], j in union(cc_lazy[1], cc_lazy[3], cc_lazy[4]) 
                                      if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)  

        con3 = JuMP.@build_constraint(sum(lom.variables[:z_var][i,j] for i in cc_lazy[3], j in union(cc_lazy[1], cc_lazy[2], cc_lazy[4]) 
                                      if !(isapprox(adjacency_augment_graph[i,j], 0, atol=1E-6))) >= 1)  
                                        
        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con1)       
        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con2)
        MOI.submit(lom.model, MOI.LazyConstraint(cb_cuts), con3)
        
        if lom.data["lazycuts_logging"]
            Memento.info(_LOGGER, "Polyhedral relaxation cuts: graph connectivity (#cc = $(length(cc_lazy)))")
        end   
         
    elseif length(cc_lazy) >= 5
        Memento.info(_LOGGER, "Polyhedral relaxation cuts: flow cuts not added for integer solutions with $(length(cc_lazy)) connected components")
        
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