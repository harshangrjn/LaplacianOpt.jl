#---------------------------------------------------------#
# Initialize all variables of the LaplacianOptModel here  #
#---------------------------------------------------------#

function variable_lifted_W_matrix(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]

    lom.variables[:W_var] = JuMP.@variable(lom.model, W_var[1:n, 1:n], Symmetric)

    for i = 1:n
        JuMP.set_lower_bound(W_var[i,i], 0)
    end

    return
end

function variable_edge_onoff(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]

    lom.variables[:z_var] = JuMP.@variable(lom.model, z_var[1:n, 1:n], Bin, Symmetric)

    return
end

function variable_algebraic_connectivity(lom::LaplacianOptModel)

    lom.variables[:γ_var] = JuMP.@variable(lom.model, γ_var >= 0)

    return
end

function variable_multi_commodity_flow(lom::LaplacianOptModel)
    n = lom.data["num_nodes"]

    lom.variables[:flow_var] = JuMP.@variable(lom.model, 0 <= flow_var[1:n,1:n,1:(n-1)] <= 1)

    return
end