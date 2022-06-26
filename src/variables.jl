#---------------------------------------------------------#
# Initialize all variables of the LaplacianOptModel here  #
#---------------------------------------------------------#

function variable_lifted_W_matrix(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]
    adjacency_base_graph = lom.data["adjacency_base_graph"]

    lom.variables[:W_var] = JuMP.@variable(lom.model, W_var[1:num_nodes, 1:num_nodes], Symmetric)

    for i = 1:num_nodes
        JuMP.set_lower_bound(W_var[i,i], 0)
    end

    return
end

function variable_edge_onoff(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]
    adjacency_base_graph = lom.data["adjacency_base_graph"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    lom.variables[:z_var] = JuMP.@variable(lom.model, z_var[1:num_nodes, 1:num_nodes], Bin, Symmetric)
    
    for i=1:num_nodes
        for j=i:num_nodes
            # No self-loop at nodes
            if i == j
                JuMP.fix(z_var[i,i], 0)
                continue
            end

            if isapprox(adjacency_augment_graph[i,j], 0, atol = 1E-6)
                # (i,j) in base graph
                if !isapprox(adjacency_base_graph[i,j], 0, atol = 1E-6)
                    JuMP.fix(z_var[i,j], 1)
                    JuMP.fix(z_var[j,i], 1)
                else
                # (i,j) neither in base not augment graph
                    JuMP.fix(z_var[i,j], 0)
                    JuMP.fix(z_var[j,i], 0)
                end
            end
        end
    end

    return
end

function variable_algebraic_connectivity(lom::LaplacianOptModel)

    lom.variables[:γ_var] = JuMP.@variable(lom.model, γ_var >= 0)

    return
end

function variable_multi_commodity_flow(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]

    lom.variables[:flow_var] = JuMP.@variable(lom.model, 0 <= flow_var[1:num_nodes, 1:num_nodes, 1:(num_nodes-1)] <= 1)

    return
end