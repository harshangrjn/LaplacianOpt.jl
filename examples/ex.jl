using LaplacianOpt
using JuMP
# using CPLEX
using GLPK
using Memento

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
    "num_nodes" => 8,
    "instance" => 1,
    "eigen_cuts_full" => true,
    "soc_linearized_cuts" => true,
    "optimizer" => "cplex"
)

#-------------------------------------------#
#      Optimization model and solution      #
#-------------------------------------------#
lom_optimizer = get_solver(params)
result = LaplacianOpt.run_LOpt_model(params, lom_optimizer)