#----------------------------------------------------------#
# Build the objective function for LaplacianOptModel here  #
#----------------------------------------------------------#

function objective_maximize_algebraic_connectivity(lom::LaplacianOptModel)
    JuMP.@objective(lom.model, Max, lom.variables[:γ_var])
end

function objective_maximize_spanning_tree_cost(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]
    edge_weights = lom.data["edge_weights"]

    JuMP.@objective(lom.model, Max, sum(edge_weights[i,j] * lom.variables[:z_var][i,j] for i=1:(n-1), j=(i+1):n))
end