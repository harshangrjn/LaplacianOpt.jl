#----------------------------------------------------------#
# Build the objective function for LaplacianOptModel here  #
#----------------------------------------------------------#

function objective_maximize_algebraic_connectivity(lom::LaplacianOptModel)
    JuMP.@objective(lom.model, Max, lom.variables[:γ_var])
end
    