import LaplacianOpt as LOpt
using JuMP
using CPLEX
using Gurobi
#using GLPK

include("solver.jl")

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
    "num_nodes" => 8,
    "instance" => 1,
    "eigen_cuts_full" => true,
    "soc_linearized_cuts" => false,
    "optimizer" => "gurobi"
)

#-------------------------------------------#
#      Optimization model and solution      #
#-------------------------------------------#
lom_optimizer = get_solver(params)
result = LOpt.run_LOpt_model(params, lom_optimizer)