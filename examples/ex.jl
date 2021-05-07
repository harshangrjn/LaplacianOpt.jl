using LaplacianOpt
using JuMP
using CPLEX
using GLPK
using Memento

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
"num_nodes" => 10,
"instance" => 1,

"optimizer" => "cplex"
)

#-------------------------------------------#
#      Optimization model and solution      #
#-------------------------------------------#
lom_optimizer = get_solver(params)
result_lopt = LaplacianOpt.run_LOpt_model(params, lom_optimizer, visualize_solution = true, visualizing_tool = "graphviz")