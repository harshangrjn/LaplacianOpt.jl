#----------------------------------------------------------#
# Build the objective function for LaplacianOptModel here  #
#----------------------------------------------------------#

function objective_maximize_algebraic_connectivity(lom::LaplacianOptModel)
    JuMP.@objective(lom.model, Max, lom.variables[:γ_var])
end

function objective_maximize_spanning_tree_cost(lom::LaplacianOptModel)
    num_nodes = lom.data["num_nodes"]
    adjacency_augment_graph = lom.data["adjacency_augment_graph"]

    JuMP.@objective(
        lom.model,
        Max,
        sum(
            adjacency_augment_graph[i, j] * lom.variables[:z_var][i, j] for
            i in 1:(num_nodes-1), j in (i+1):num_nodes if
            adjacency_augment_graph[i, j] > 0
        )
    )
end
