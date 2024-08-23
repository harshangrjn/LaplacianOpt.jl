using LaplacianOpt
using JuMP
using Gurobi
using JSON
using LinearAlgebra
using Clarabel

include("optimizers.jl")

function constraint_lazycallback(cb_cuts)
    status = JuMP.callback_node_status(cb_cuts, model)

    if status in [MOI.CALLBACK_NODE_STATUS_INTEGER]

        X_val = JuMP.callback_value.(Ref(cb_cuts), X_psd)
        violated_eigen_vec = LaplacianOpt._violated_eigen_vector(X_val, tol = 1E-10)
        
        # Eigen-cuts in X_psd-space
        if violated_eigen_vec !== nothing
            con = JuMP.@build_constraint(violated_eigen_vec' * X_psd * violated_eigen_vec >= 0)
            MOI.submit(model, MOI.LazyConstraint(cb_cuts), con)
            # println("eigen cut added")
        end
    end
end

function constraints_boolean_quadric(model, X_psd::LinearAlgebra.Symmetric{VariableRef, Matrix{VariableRef}})
    num_nodes = size(X_psd)[1]

    JuMP.@constraint(model, [i=1:num_nodes, j=1:num_nodes], X_psd[i,j] <= X_psd[i,i])

    for i=1:num_nodes, j=1:num_nodes, l=1:num_nodes
        JuMP.@constraint(model, X_psd[i,l] + X_psd[j,l] - X_psd[i,j] <= X_psd[l,l])
        JuMP.@constraint(model, X_psd[i,i] + X_psd[j,j] - X_psd[i,j] <= 1)
        JuMP.@constraint(model, X_psd[i,i] + X_psd[j,j] + X_psd[l,l] - X_psd[i,j] - X_psd[i,l] - X_psd[j,l] <= 1)
    end

    return
end

num_nodes = 10
instance = 10
optimizer = "clarabel"
formulation = "set_cardinality" # set_volume, set_cardinality
boolean_quadric = true # additional facets

data_dict = LaplacianOpt.parse_file("instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
sol_dict = JSON.parsefile("solutions/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")

adjacency_matrix = zeros(num_nodes, num_nodes)

for k = 1:length(sol_dict["graph_adjacency"])
    I = sol_dict["graph_adjacency"][k][1]
    J = sol_dict["graph_adjacency"][k][2]
    adjacency_matrix[I,J] = data_dict["adjacency_augment_graph"][I,J]
    adjacency_matrix[J,I] = adjacency_matrix[I,J]
end

if optimizer == "clarabel"
    model = JuMP.Model(Clarabel.Optimizer)
    JuMP.set_optimizer_attribute(model, "verbose", true)
    JuMP.set_optimizer_attribute(model, "equilibrate_enable",false)
else
    model = JuMP.Model(get_gurobi(solver_log = true))
end

# Variables
JuMP.@variable(model, 0 <= X_psd[1:num_nodes, 1:num_nodes] <= 1, Symmetric)
JuMP.@variable(model, phi >= 0)
if optimizer !== "clarabel"
    JuMP.@variable(model, 0 <= z_dummy <= 1, Bin)
end

# Constraints
expr = AffExpr()
for i=1:(num_nodes-1), j=(i+1):num_nodes
    if !isapprox(adjacency_matrix[i,j], 0, atol = 1E-6)
        global expr += adjacency_matrix[i,j] * (X_psd[i,i] + X_psd[j,j] - 2*X_psd[i,j])
    end
end

JuMP.@constraint(model, phi >= expr)

if formulation == "set_cardinality"
    JuMP.@constraint(model, sum(X_psd) >= 1)
    JuMP.@constraint(model, sum(X_psd) <= num_nodes/2)
    JuMP.@constraint(model, sum(X_psd[i,i] for i=1:num_nodes) == 1)

elseif formulation == "set_volume"

    # Laplacian to get unweighted volume (sum of degrees of nodes) of the set
    laplacian = LaplacianOpt.laplacian_matrix(float(adjacency_matrix .> 0))
    
    # k_lb = minimum(unique(float(adjacency_matrix .> 0))[2:end])
    k_lb = 1
    k_ub = sum(adjacency_matrix .> 0) / 2 # |E|

    JuMP.@constraint(model, sum(X_psd[i,i]*laplacian[i,i] for i=1:num_nodes) == 1)

    expr2 = AffExpr()
    for i = 1:num_nodes
        global expr2 += laplacian[i,i]^2 * X_psd[i,i]
        if i <= num_nodes - 1
            for j = (i+1):num_nodes
                expr2 += 2 * X_psd[i,j] * laplacian[i,i] * laplacian[j,j]
            end
        end
    end

    JuMP.@constraint(model, expr2 >= k_lb)
    JuMP.@constraint(model, expr2 <= k_ub)
end

if optimizer == "clarabel"
    JuMP.@constraint(model, X_psd >= 0, PSDCone())
else
    # Callback for eigen cuts
    JuMP.set_attribute(model, MOI.LazyConstraintCallback(), constraint_lazycallback)    
end

# Boolean quadric polytope facets
boolean_quadric && constraints_boolean_quadric(model, X_psd)

# Objective
JuMP.@objective(model, Min, phi)

JuMP.optimize!(model)

println("=====================================")
println("Objective value: ", JuMP.objective_value(model))
println("Lambda_2/2: ", sol_dict["lambda_2"]/2)
println("Cheeger constant (set cardinality based): ", LaplacianOpt.cheeger_constant(adjacency_matrix, get_gurobi())["cheeger_constant"])
println("=====================================")

