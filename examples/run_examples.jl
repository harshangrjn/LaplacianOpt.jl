import LaplacianOpt as LOpt
using JuMP
using Gurobi
# using GLPK
# using CPLEX

include("optimizer.jl")

lopt_optimizer = get_gurobi()

num_nodes = 8
instance  = 1

# Data format has to be as given in this JSON file
file_path = joinpath(@__DIR__, "instances/$(num_nodes)_nodes/$(num_nodes)_$(instance).json")
data_dict = LOpt.parse_file(file_path)

#-------------------------------#
#      User-defined inputs      #
#-------------------------------#
params = Dict{String, Any}(
    "data_dict"           => data_dict,
    "augment_budget"      => (num_nodes-1),
    "eigen_cuts_full"     => true,
    "soc_linearized_cuts" => false,
    "eigen_cuts_2minors"  => false,
    "eigen_cuts_3minors"  => false
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