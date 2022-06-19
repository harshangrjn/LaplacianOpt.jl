import LaplacianOpt as LOpt
using JuMP
# using CPLEX
using Gurobi
#using GLPK

include("optimizer.jl")

result_all = zeros(10,3)
lopt_optimizer = get_gurobi()

num_nodes = 8
instance = 1

# Data format has to be as given in this JSON file
file_path = joinpath(dirname(pathof(LOpt)),"..", "examples/instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
data_dict = LOpt.parse_file(file_path)

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
    "instance"            => instance,
    "data_dict"           => data_dict,
    "augment_budget"      => (num_nodes-1),
    "eigen_cuts_full"     => true,
    "soc_linearized_cuts" => false
    )

#-------------------------------------------#
#      Optimization model and solution      #
#-------------------------------------------#
result = LOpt.run_LOpt(params, lopt_optimizer)

#------------------------------------#
#      Visualize solution graph      #
#------------------------------------#
# data = LOpt.get_data(params)
# LOpt.visualize_solution(result, data; display_edge_weights = false)
