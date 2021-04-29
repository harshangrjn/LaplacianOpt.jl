#------------------------------------------------------#
# Build all constraints of the LaplacianOptModel here  #
#------------------------------------------------------#

function constraint_build_W_var_matrix(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]
    ew_matrix = lom.data["edge_weights"]

    # Diagonal entries
    JuMP.@constraint(lom.model, [i=1:n], lom.variables[:W_var][i,i] == sum(ew_matrix[i,:] .* lom.variables[:z_var][i,:]) - lom.variables[:γ_var]*(n-1)/n)

    # Off-diagonal entries
    JuMP.@constraint(lom.model, [i=1:(n-1), j=(i+1):n], lom.variables[:W_var][i,j] == -ew_matrix[i,j] * lom.variables[:z_var][i,j] + lom.variables[:γ_var]/n)

    return
end

function constraint_topology_no_self_loops(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]
    
    JuMP.@constraint(lom.model, [i=1:n], lom.variables[:z_var][i,i] == 0)
    
    return
end

function constraint_topology_vertex_cutset(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]
    
    JuMP.@constraint(lom.model, [i=1:n], sum(lom.variables[:z_var][i,j] for j in setdiff((1:n), i)) >= 1)
    
    return
end

function constraint_topology_total_edges(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]

    JuMP.@constraint(lom.model, sum(lom.variables[:z_var][i,j] for i=1:(n-1), j=(i+1):n) == (n-1))
    
    return
end